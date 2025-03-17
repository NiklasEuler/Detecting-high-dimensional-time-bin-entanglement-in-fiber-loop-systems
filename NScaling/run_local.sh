#!/bin/bash
cd "$(dirname "$0")"

export RUNNAME=RUN_LOCAL

mkdir -p ${RUNNAME}/Logs
mkdir -p ${RUNNAME}/Data

cp run_local.sh ${RUNNAME}/
cp run.jl ${RUNNAME}/

cd ${RUNNAME}/

export N_RESAMPLES=500 # number of independent experimental runs
export N_MAX=16 # maximum number of time bins
export EPS_ANGL=0.02 # angular deviation on the central coupler in units of pi

export N_SAMPLES_POP="[500,950,1950,3975]" # number of samples for the population measurement for each N
export N_SAMPLES_PHASE_EACH="[500,1000,1975,3950]" # number of samples for each phase estimation DTQW for each N
export N_SAMPLES_COMPOUND_EACH="[500,350,300,275]" # number of samples for each of the compound-setting DTQWs for each N

export T=4 # number of Julia threads
export USE_MKL=false # flag to use MKL or OPENBLAS
export MKL_T=1 # number of MKL threads
export BLAS_T=1 # number of OPENBLAS threads
julia -t ${T} run.jl $USE_MKL $N_RESAMPLES $N_MAX $EPS_ANGL $N_SAMPLES_POP $N_SAMPLES_PHASE_EACH $N_SAMPLES_COMPOUND_EACH 2>&1 | tee Logs/Log.txt
