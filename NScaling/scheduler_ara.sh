#!/bin/bash

export RUNNAME=Rerun

mkdir -p ${RUNNAME}/Logs
mkdir -p ${RUNNAME}/Data

cp startup_ara.sbatch ${RUNNAME}/
cp run.jl ${RUNNAME}/
cp scheduler_ara.sh ${RUNNAME}/

cd ${RUNNAME}/

export N_RESAMPLES=500 # number of independent experimental runs
export N_MAX=16 # maximum number of time bins
export EPS_ANGL=0.02 # angular deviation on the central coupler in units of pi

export N_SAMPLES_POP="[500,950,1950,3975]" # number of samples for the population measurement for each N
export N_SAMPLES_PHASE_EACH="[500,1000,1975,3950]" # number of samples for each phase estimation DTQW for each N
export N_SAMPLES_COMPOUND_EACH="[500,350,300,275]" # number of samples for each of the compound-setting DTQWs for each N

export T=72
export CORES="--cpus-per-task=$T"
export MEM="--mem=96gb"
export TIME="--time=03:00:00"
export MKL_T=1
export BLAS_T=1
sbatch --job-name=${RUNNAME}_T${T}LT${MKL_T}MKL ${MEM} ${TIME} ${CORES} --output=Logs/log_out_T${T}LT${MKL_T}MKL.txt startup_ara.sbatch $T $MKL_T false $N_RESAMPLES $N_MAX $EPS_ANGL $N_SAMPLES_POP $N_SAMPLES_PHASE_EACH $N_SAMPLES_COMPOUND_EACH
sleep 0.1



