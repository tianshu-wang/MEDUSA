#!/bin/bash
@(#SBATCH -d afterok:@CHAINED_JOB_ID@)@
#SBATCH -p regular
#SBATCH -t @WALLTIME@
#SBATCH -N @NODES@ -n @NUM_PROCS@ -c @NCPU_PER_PROC@
#SBATCH -C knl,quad,cache
#SBATCH -L SCRATCH
#SBATCH -J @SIMULATION_NAME@
#SBATCH -e @RUNDIR@/@SIMULATION_NAME@.err
#SBATCH -o @RUNDIR@/@SIMULATION_NAME@.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=@EMAIL@

module load cray-hdf5-parallel

echo "Preparing:"
set -x                          # Output commands
set -e                          # Abort on errors

cd @RUNDIR@

echo "Checking:"
pwd
hostname
date

echo "Environment:"
export GMON_OUT_PREFIX=gmon.out
export OMP_NUM_THREADS=@NUM_THREADS@
export OMP_PROC_BIND=true
export OMP_PLACES=threads
env | sort > ENVIRONMENT
echo ${SLURM_NODELIST} > NODES

sbcast --compress=lz4 @EXECUTABLE@ /tmp/fornax.exe

echo "Starting:"
time srun \
    -N @NODES@ \
    -n @NUM_PROCS@ \
    -c @NCPU_PER_PROC@ \
    --cpu_bind=@CPU_BIND@ \
    /tmp/fornax.exe \
    @PARFILE@

echo "Stopping:"
date

echo "Done."
