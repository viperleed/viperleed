#!/bin/bash

#SBATCH -J viperleed_testv1
# use 5 cores
#SBATCH -n 5
#SBATCH --mem=8G

#SBATCH --qos=devel_0096
#SBATCH --partition=mem_0096
#SBATCH --account=p71704

#SBATCH --ntasks-per-core=2


# load conda environment - must already exist
module purge
spack load -r /fffbmf3
conda activate viperleed

# intel compilers
module load intel/19.1.3 intel-mpi/2019.7

echo "Loading finished"

vpr_path=$DATA/source/
work_path=$DATA/Tests/

echo "ViPErLEED Source:		$vpr_path"
echo "Working directory:	$work_path"
echo

# run ViPErLEED in Singularity image – /gpfs needs to be bound to /gpfs for access to $DATA
cd $work_path
singularity exec --bind /gpfs:/gpfs $image_path python3 job.py -s $vpr_path -w $work_path

echo "ViPErLEED finished"
