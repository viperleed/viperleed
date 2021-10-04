# -*- coding: utf-8 -*-
"""
@author: Alexander M. Imre

Use this to recompile some precompiled executables when moving to a new workstation.
In particular, EEASiSSS, beamgen and random_.c need to be compiled.
"""

import subprocess
from pathlib import Path
import os
import logging

############################
# SPECIFY COMPILERS						# OpenMPI or Intel compilers
############################
TensErLEED_version = "1.71"

c_compiler = "gcc"						# "gcc" or "icc"
fortran_compiler = "gfortran"			# "gfortran" or "ifort"

c_mpi_compiler = "mpicc"				# "mpicc" or "mpiicc"
fortran_mpi_compiler = "mpifort"		# "mpifort" or "mpiifort"

suppress_fortran_warnings = True		# Supress warnings in log file. (Deprecated features in soure will throw a large amount of warnings and make logfile hard to read.)

############################
# SCRIPT BELOW
############################

class to_compile:
	logname = "" #class variable; shared by all instances!
	
	def __init__(self):
		self.compiler = ""
		self.sources = []
		self.destination = ""
		self.flags = ""
		self.logname = ""
	
	def set_compiler(self, compiler):
		self.compiler = compiler
		
	def set_flags(self, flags):
		self.flags = flags
		
	def set_logname(self, logname):
		self.logname = logname
		
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
		if os.path.isfile(self.logname):
			sep = "\n\n"
		try:
			with open(self.logname, "a") as log:
				log.write(sep + "############\n# COMPILING: " + sources_str
							+ "\n############\n\n")
			with open(self.logname, "a") as log:
				r = subprocess.run(compile_command.split(), stdout=log, stderr=log)
				log.write("Command sent:\n" + compile_command + "\n")
				if self.flags == "":
					log.write(sep + "No compile flags were set. Intentional?")
		except Exception:
			logger.error("Error compiling "+ sources_str)
			raise
		if r.returncode != 0:
			raise RuntimeError("Compiler subprocess returned {}"
								.format(r.returncode))
		
		return None
	
	
###############################################
#                  MAIN                       #
###############################################
def main():
	
	logname="recompile_parts.log" # set logfile
	logger = logging.getLogger('recompile_parts')
	
	this_directory = Path(".").resolve()

	############################
	# random
	############################
	
	source_filename = "random_.c"
	destination_filename = "random_.o"
	flags = "-c -Ofast"
	
	source_path = this_directory.parent / 'tensorleed' / ('TensErLEED-v' + TensErLEED_version) / 'lib' / source_filename
	destination_path = this_directory.parent / 'tensorleed' / ('TensErLEED-v' + TensErLEED_version) / 'lib' / destination_filename
	sources = [source_path]
	
	# is needed with and without MPI (TODO: why?)
	random = to_compile()
	random.set_compiler(c_compiler)
	random.set_logname(logname)
	random.set_destination(destination_path)
	random.add_source(source_path)
	random.set_flags(flags)
	random.send_compile(suppress_fortran_warnings) # compile without MPI
	
	# below: compile with MPI
	random.set_compiler(c_mpi_compiler)
	destination_filename = "MPIrandom_.o"
	destination_path = this_directory.parent / 'tensorleed' / ('TensErLEED-v' + TensErLEED_version) / 'lib' / destination_filename
	random.set_destination(destination_path)
	random.send_compile(suppress_fortran_warnings)
	
	
	############################
	# EEASiSSS
	############################
	
	EEASiSSS = to_compile()
	sources = ["eeasisss.f90", "imported_routines.f90"]
	destination_filename = "EEASiSSS.x"
	flags = "-Ofast"									# optimization level 3
	EEASiSSS.set_compiler(fortran_mpi_compiler) 		# This needs to be MPI compiled
	
	for source_filename in sources:
		source_path = this_directory.parent / 'tensorleed' / 'eeasisss_code' / 'modified' / source_filename
		EEASiSSS.add_source(source_path)
	destination_path = this_directory.parent / 'tensorleed' / destination_filename
	EEASiSSS.set_destination(destination_path)
	EEASiSSS.set_logname(logname)
	
	# Place module (.mod) files in specified directory rather then working directory (= utilities directory)
	# See here for details: https://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules
	if EEASiSSS.compiler == 'mpifort':
		flags += " -J" + str((this_directory.parent / 'tensorleed' / 'eeasisss_code' / 'modified' / 'mod_files'))
	elif EEASiSSS.compiler == 'mpiifort':
		flags += " -module " + str((this_directory.parent / 'tensorleed' / 'eeasisss_code' / 'modified' / 'mod_files'))
	EEASiSSS.set_flags(flags)
	EEASiSSS.send_compile(suppress_fortran_warnings)
	
	############################
	# beamgen
	############################
	
	beamgen_version = "1.7"
	
	beamgen = to_compile()
	sources = ["beamgen.v" +beamgen_version + ".f"]
	destination_filename = "beamgen.v" +beamgen_version 	# no file ending!
	flags = "-Ofast"										# optimization level 3
	beamgen.set_compiler(fortran_mpi_compiler)				# This needs to be MPI compiled
	
	for source_filename in sources:
		source_path = this_directory.parent / 'tensorleed' / 'beamgen_source' / source_filename
		beamgen.add_source(source_path)
	destination_path = this_directory.parent / 'tensorleed' / destination_filename
	beamgen.set_destination(destination_path)
	beamgen.set_logname(logname)
	
	beamgen.set_flags(flags)
	beamgen.send_compile(suppress_fortran_warnings)

	############################
	# EEASISSS -new
	############################

	EEASiSSS = to_compile()
	sources = ["eeasisss_init.f90", "eeasisss_main.f90", "eeasisss.f90"]
	destination_filename = "eeasisss"
	flags = "-Ofast"  # optimization level 3
	EEASiSSS.set_compiler(fortran_mpi_compiler)  # This needs to be MPI compiled

	for source_filename in sources:
		source_path = this_directory.parent / 'tensorleed' / 'eeasisss_new' / source_filename
		EEASiSSS.add_source(source_path)
	destination_path = this_directory.parent / 'tensorleed' / 'eeasisss_new' /destination_filename
	EEASiSSS.set_destination(destination_path)
	EEASiSSS.set_logname(logname)

	# Place module (.mod) files in specified directory rather then working directory (= utilities directory)
	# See here for details: https://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules
	if EEASiSSS.compiler == 'mpiifort':
		flags += " -I " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
		flags += " -module" + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
	elif EEASiSSS.compiler == 'mpifort':
		flags += " -I " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
		flags += " -J " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
	EEASiSSS.set_flags(flags)
	EEASiSSS.send_compile(suppress_fortran_warnings)

	EEAS = to_compile()
	sources = ["eeas.f90"]
	destination_filename = "eeas"
	flags = "-Ofast"  # optimization level 3
	EEAS.set_compiler(fortran_compiler)  # This needs to be MPI compiled

	for source_filename in sources:
		source_path = this_directory.parent / 'tensorleed' / 'eeasisss_new' / source_filename
		EEAS.add_source(source_path)
	destination_path = this_directory.parent / 'tensorleed' / 'eeasisss_new' / destination_filename
	EEAS.set_destination(destination_path)
	EEAS.set_logname(logname)

	# Place module (.mod) files in specified directory rather then working directory (= utilities directory)
	# See here for details: https://stackoverflow.com/questions/8855896/specify-directory-where-gfortran-should-look-for-modules
	if EEAS.compiler == 'ifort':
		flags += " -I " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
		flags += " -module" + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
	elif EEAS.compiler == 'gfortran':
		flags += " -I " + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
		flags += " -J" + str((this_directory.parent / 'tensorleed' / 'eeasisss_new' / 'mod_files'))
	EEAS.set_flags(flags)
	EEAS.send_compile(suppress_fortran_warnings)

if __name__ == "__main__":
		main()