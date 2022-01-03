#!/bin/bash

#SBATCH -J viperleed_testv1
# use 5 cores
#SBATCH -n 5
#SBATCH --mem=8G

#SBATCH --qos=mem_0096 # can put to devel_0096 for testing
#SBATCH --partition=mem_0096
#SBATCH --account=p71704

#SBATCH --ntasks-per-core=2


# load conda environment - must already exist
module purge
spack load -r /fffbmf3 # miniconda3
source ~/.bashrc
conda activate viperleed

# intel compilers - DO NOT USE intel-mpi/2019 # A known issue will cause error!
module load intel/19.1.3 intel-mpi/latest
module load intel-mkl/2019.3
echo "Loading finished"

vpr_path=$DATA/source/
work_path=$DATA/Tests/

echo "ViPErLEED Source:         $vpr_path"
echo "Working directory:        $work_path"
echo

# run ViPErLEED in conda environment with all required packages
cd $work_path
python3 job.py -s $vpr_path -w $work_path

echo "ViPErLEED finished"
