# -*- coding: utf-8 -*-
"""
@author: Alexander M. Imre

Use this to recompile some precompiled executables when moving to a new workstation or changes to the environment are made.
In particular, EEASiSSS, beamgen and random_.c need to be compiled.
"""

import subprocess
from pathlib import Path
import os
import logging
import argparse
import time

suppress_fortran_warnings = True		# Supress warnings in log file. (Deprecated features in soure will throw a large amount of warnings and make logfile hard to read.)

############################
# SCRIPT BELOW
############################

class to_compile:
	
	def __init__(self, logger):
		self.compiler = ""
		self.sources = []
		self.destination = ""
		self.flags = ""
		self.logger = logger
	
	def set_compiler(self, compiler):
		self.compiler = compiler
		
	def set_flags(self, flags):
		self.flags = flags

		
	def add_source(self, source_file):
		self.sources.append(source_file)
		
	def set_destination(self, destination_file):
		self.destination = destination_file
		
	def send_compile(self, suppress_fortran_warnings):
		# adapted from leedbase.py by Florian Kraushofer
		sources_str = ' '.join([str(i) for i in self.sources])
		if suppress_fortran_warnings:
			if (self.compiler == "gfortran") or (self.compiler == "mpifort"):
				self.flags += " -w"
			elif (self.compiler == "ifort") or (self.compiler == "mpiifort"):
				self.falgs += " -nowarn"
				
		compile_command = self.compiler + " -o " + str(self.destination) + ' ' + sources_str +  ' ' + self.flags
		
		sep = ""
		try:
			self.logger.info(sep + "COMPILING: " + sources_str)
			self.logger.info("Command sent:\n" + compile_command + "\n")
			r = subprocess.run(compile_command.split(), capture_output=True, text=True)
			if self.flags == "":
				self.logger.warning(sep + "No compile flags were set. Intentional?")
			if r.stdout != "":
				self.logger.info(r.stdout)
			if r.stderr != "":
				self.logger.error(r.stderr)
		except Exception:
			self.logger.error("Error compiling " + sources_str)
			raise
		if r.returncode != 0:
			raise RuntimeError(f"Compiler subprocess returned {r.returncode}")
		else:
			self.logger.info("Compiler returncode 0")

		
		return None
	
	
###############################################
#                  MAIN                       #
###############################################
def main():

	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-s", "--source",
		help=("specify ViPErLEED source directory (without final /viperleed)"),
		type=str)
	parser.add_argument(
		"-l", "--logfile",
		help=("specify logfile location"),
		type=str)
	parser.add_argument(
		"-c", "--compiler",
		help=("specify which compiler to use"),
		type=str)
	parser.add_argument(
		"-v", "--version",
		help=("specify for which TensErLEED version to compile"),
		type=str)

	args = parser.parse_args()

	this_directory = Path(".").resolve()

	if args.source:
		source_directory = args.source
		# check if it exists
		if not (os.path.exists(source_directory) and os.path.exists(os.path.join(source_directory, "viperleed"))):
			raise Exception(f"Directory {source_directory} not found")
	else:
		# no source directory specified
		raise Exception("No source directory specified, use -s option")

	if args.logfile:
		logname = args.logfile
	else:
		timestamp = time.strftime("%y%m%d-%H%M%S", time.localtime())
		logname = 'recompile_parts-' + timestamp + '.log'

	logger = logging.getLogger("recompile_parts")
	# start logger, write to file:
	logger.setLevel(logging.INFO)
	logFormatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fileHandler = logging.FileHandler(logname, mode="w")
	fileHandler.setFormatter(logFormatter)
	logger.addHandler(fileHandler)

	logger.info("Recompile script started.")

	# Set compilers
	if ((args.compiler == "gcc") or (args.compiler == "gnu") or (args.compiler == "g")):
		comp = "gnu"
		c_compiler = "gcc"
		fortran_compiler = "gfortran"
		c_mpi_compiler = "mpicc"
		fortran_mpi_compiler = "mpifort"
	elif ((args.compiler == "intel") or (args.compiler == "i") or (args.compiler == "icc")):
		comp = "intel"
		c_compiler = "icc"
		fortran_compiler = "ifort"
		c_mpi_compiler = "mpiicc"
		fortran_mpi_compiler = "mpiifort"
	else:
		raise Exception("No compiler specified, use -c option")

	# Set version to use
	if args.version:
		# check if valid version
		TensErLEED_directoy = os.path.join(source_directory, "viperleed", 'tensorleed', 'TensErLEED-v' + args.version)
		if not os.path.exists(TensErLEED_directoy):
			raise Exception(f"Invalid TensErLEED version, not directory {TensErLEED_directoy}")
		TensErLEED_path = Path(TensErLEED_directoy)
	else:
		raise Exception("No TensErLEED version specified, use -v option")


	############################
	# random
	############################
	
	source_filename = "random_.c"
	destination_filename = "random_.o"
	flags = "-c -O2"
	
	source_path = TensErLEED_path / 'lib' / source_filename
	destination_path = TensErLEED_path / 'lib' / destination_filename
	sources = [source_path]
	
	# is needed with and without MPI (TODO: why?)
	random = to_compile(logger)
	random.set_compiler(c_compiler)
	random.set_destination(destination_path)
	random.add_source(source_path)
	random.set_flags(flags)
	random.send_compile(suppress_fortran_warnings) # compile without MPI
	
	# below: compile with MPI
	random.set_compiler(c_mpi_compiler)
	destination_filename = "MPIrandom_.o"
	destination_path = TensErLEED_path / 'lib' / destination_filename
	random.set_destination(destination_path)
	random.send_compile(suppress_fortran_warnings)
	
	
	############################
	# EEASiSSS
	############################
	
	EEASiSSS = to_compile(logger)
	sources = ["eeasisss.f90", "imported_routines.f90"]
	destination_filename = "EEASiSSS.x"
	flags = "-O2"
	EEASiSSS.set_compiler(fortran_compiler) 		# This needs to be MPI compiled
	
	for source_filename in sources:
		source_path = this_directory.parent / 'tensorleed' / 'eeasisss_code' / 'modified' / source_filename
		EEASiSSS.add_source(source_path)
	destination_path = this_directory.parent / 'tensorleed' / destination_filename
	EEASiSSS.set_destination(destination_path)
	
	# Place module (.mod) files in specified directory rather then working directory (= utilities directory)
	# See here for details: https://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules
	if comp == "gnu":
		flags += " -J " + str((this_directory.parent / 'tensorleed' / 'eeasisss_code' / 'modified' / 'mod_files'))
	elif comp == "intel":
		flags += " -module " + str((this_directory.parent / 'tensorleed' / 'eeasisss_code' / 'modified' / 'mod_files'))
	else:
		raise Exception("Sanity check at EEASISSS compiler selection failed")
	EEASiSSS.set_flags(flags)
	EEASiSSS.send_compile(suppress_fortran_warnings)
	
	############################
	# beamgen
	############################
	
	beamgen_version = "1.7"
	
	beamgen = to_compile(logger)
	sources = ["beamgen.v" +beamgen_version + ".f"]
	destination_filename = "beamgen.v" +beamgen_version 	# no file ending!
	flags = "-O2"										# optimization level 3
	beamgen.set_compiler(fortran_mpi_compiler)				# This needs to be MPI compiled
	
	for source_filename in sources:
		source_path = this_directory.parent / 'tensorleed' / 'beamgen_source' / source_filename
		beamgen.add_source(source_path)
	destination_path = this_directory.parent / 'tensorleed' / destination_filename
	beamgen.set_destination(destination_path)

	beamgen.set_flags(flags)
	beamgen.send_compile(suppress_fortran_warnings)

"""
Unused for now!
	############################
	# EEASISSS -new
	############################

	EEASiSSS = to_compile()
	sources = ["eeas.f90", "eeasisss_init.f90", "eeasisss_main.f90", "eeasisss.f90"] # Order important since they depended on modules in previous files!!
	destination_filename = "eeasisss"
	flags = "-O2"  # safe optimization
	EEASiSSS.set_compiler(fortran_compiler)

	for source_filename in sources:
		source_path = this_directory.parent / 'tensorleed' / 'eeasisss_new' / source_filename
		EEASiSSS.add_source(source_path)
	destination_path = this_directory.parent / 'tensorleed' / 'eeasisss_new' /destination_filename
	EEASiSSS.set_destination(destination_path)

	# Place module (.mod) files in specified directory rather then working directory (= utilities directory)
	# See here for details: https://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules
	if EEASiSSS.compiler == 'ifort':
		flags += " -I " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
		flags += " -module " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
	elif EEASiSSS.compiler == 'gfortran':
		flags += " -I " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
		flags += " -J " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
	EEASiSSS.set_flags(flags)
	EEASiSSS.send_compile(suppress_fortran_warnings)
"""

if __name__ == "__main__":
		main()