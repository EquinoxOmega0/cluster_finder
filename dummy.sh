#!/bin/sh
#PBS -N make_mock_catalogue
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1000:00:00
cd $PBS_O_WORKDIR
export PATH=$PATH:$PBS_O_WORKDIR
export LD_LIBRARY_PATH=/your_libraries_path:$LD_LIBRARY_PATH

REQUIRED_CPUS="`cat $PBS_NODEFILE | wc -l`"

mpiexec -n ${REQUIRED_CPUS} ./a.out
