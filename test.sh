### Job Name
#PBS -N PARALLEL_GRAPE_JOB
### Submits the job to the default queue
#PBS -q workq
### Requests 4 CPUs on different nodes
#PBS -l nodes=n01.local:ppn=1+n02.local:ppn=1+n03.local:ppn=1+n04.local:ppn=1
### Output files
#PBS -o par_grape_pbs.log
#PBS -e par_grape_pbs.err
### Declare job non-rerunable
#PBS -r n

#!/bin/sh
# Change to working directory
cd $PBS_O_WORKDIR
# Setting environment variable
export PATH=:$PBS_O_WORKDIR

echo "Job started at `date`"
echo ""
# Alternative variant to submit grape parallel jobs
# Run your executable program (full path or your name of pogram in working directory)
# MVAPICH2 case
mpiexec-mvapich2 -n 4 ./a.out
echo ""
echo "Job ended at `date`"