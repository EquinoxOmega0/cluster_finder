#PBS -N mcc
#PBS -q workq
#PBS -l nodes=1:ppn=8
#PBS -o test_pbs.log
#PBS -e test_pbs.err
#PBS -r n

#!/bin/sh

# Change to working directory
cd $PBS_O_WORKDIR
# Setting environment variable
export PATH=:$PBS_O_WORKDIR

echo "Job started at `date`"
echo ""
# Run your executable program (full path or your name of program in working directory)
mpiexec-mvapich2 -n 1 ./a.out
echo ""
echo "Job ended at `date`"