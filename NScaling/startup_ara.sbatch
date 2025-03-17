#!/bin/bash

#SBATCH --partition=s_test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

#SBATCH --cpus-per-task=72

module load julia
julia -v
pwd
export MKL_DYNAMIC=false # Disable dynamic MKL thread number adjustment
export JULIA_NUM_THREADS=$1 # Number of Julia threads
export OPENBLAS_NUM_THREADS=$2 # Number of OpenBLAS threads
export MKL_NUM_THREADS=$2 # Number of MKL threads

export USE_MKL=$3 # false: OpenBLAS, true: MKL
export N_RESAMPLES=$4 # Number of resamples
export N_MAX=$5 # Maximum number of time bins
export EPS_ANGL=$6 # angular deviation on the central coupler in units of pi

export N_SAMPLES_POP=$7 # Array of number of samples for the population measurement for all N
export N_SAMPLES_PHASE_EACH=$8 # Array of number of samples for each phase measurement DTQW for all N
export N_SAMPLES_COMPOUND_EACH=$9 # Array of number of samples for each compound-setting DTQW for all N


srun julia ./run.jl $USE_MKL $N_RESAMPLES $N_MAX $EPS_ANGL $N_SAMPLES_POP $N_SAMPLES_PHASE_EACH $N_SAMPLES_COMPOUND_EACH


