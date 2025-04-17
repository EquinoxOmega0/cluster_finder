#!/bin/sh
#PBS -N mill_vis
#PBS -q new
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
cd $PBS_O_WORKDIR
export PATH=$PATH:$PBS_O_WORKDIR
export LD_LIBRARY_PATH=/your_libraries_path:$LD_LIBRARY_PATH

REQUIRED_CPUS="`cat $PBS_NODEFILE | wc -l`"

mpirun -n ${REQUIRED_CPUS} ./a.out
