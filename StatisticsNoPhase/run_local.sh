#!/bin/bash
cd "$(dirname "$0")"

export RUNNAME=RUN_LOCAL

mkdir -p ${RUNNAME}/Logs
mkdir -p ${RUNNAME}/Data

cp run_local.sh ${RUNNAME}/
cp run.jl ${RUNNAME}/

cd ${RUNNAME}/

export N_RESAMPLES=500 # number of independent experimental runs
export N=8 # number of time bins
export N_SAMPLES=8000 # total number of samples
export EPS_ANGL=0.02 # angular deviation on the central coupler in units of pi
export EPS=0.05215 # dephasing parameter
# 0.06905 for N = 2 | 0.054836 for N = 4 | 0.05215 for N = 8 | 0.05152 for N = 16

export N_SAMPLES_COMPOUND_EACH="[140,200,260,320,380,440,500,560]"
    # number of samples for each of the compound-setting DTQWs

export JULIA_NUM_THREADS=4 # number of Julia threads
export USE_MKL=false # flag to use MKL or OPENBLAS
export MKL_T=1 # number MKL threads
export BLAS_T=1 # number of OPENBLAS threads
julia run.jl $USE_MKL $N_RESAMPLES $N $N_SAMPLES $EPS_ANGL $EPS $N_SAMPLES_COMPOUND_EACH 2>&1 | tee Logs/Log.txt



