#!/bin/bash

#SBATCH --partition=s_standard
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

module load julia
julia -v
pwd
export MKL_DYNAMIC=false # Disable dynamic MKL thread number adjustment
export JULIA_NUM_THREADS=$1 # Number of Julia threads
export OPENBLAS_NUM_THREADS=$2 # Number of OpenBLAS threads
export MKL_NUM_THREADS=$2 # Number of MKL threads

export USE_MKL=$3 # false: OpenBLAS, true: MKL
export N_RESAMPLES=$4 # Number of resamples
export N=$5 # Number of time bins
export N_SAMPLES=$6 # Number of total samples in each experimental run
export EPS_ANGL=$7 # angular deviation on the central coupler in units of pi
export EPS=$8 # dephasing parameter

export N_SAMPLES_COMPOUND_EACH=$9 # Array of number of samples for each compound-setting DTQW

srun julia ./run.jl $USE_MKL $N_RESAMPLES $N $N_SAMPLES $EPS_ANGL $EPS $N_SAMPLES_COMPOUND_EACH


