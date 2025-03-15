#!/bin/bash

export RUNNAME=Rerun

mkdir -p ${RUNNAME}/Data
mkdir -p ${RUNNAME}/Logs
mkdir -p ${RUNNAME}/Figs

cp startup_ara.sbatch ${RUNNAME}/
cp run.jl ${RUNNAME}/
cp scheduler_ara.sh ${RUNNAME}/

cd ${RUNNAME}/

export N_RESAMPLES=500 # number of independent experimental runs
export N=8 # number of time bins
export N_SAMPLES=4000 # total number of samples
export EPS_ANGL=0.02 # angular deviation on the central coupler in units of pi
export EPS=0.05215 # dephasing parameter
# 0.06905 for N = 2 | 0.054836 for N = 4 | 0.05215 for N = 8 | 0.05152 for N = 16

export N_SAMPLES_COMPOUND_EACH="[140,200,260,320,380,440,500,560]"
    # number of samples for each of the compound-setting DTQWs

export T=72
export CORES="--cpus-per-task=$T"
export MEM="--mem=96gb"
export TIME="--time=03:00:00"
export MKL_T=1
export BLAS_T=1
sbatch --job-name=${RUNNAME}_T${T}LT${MKL_T}MKL ${MEM} ${TIME} ${CORES} --output=Logs/log_out_T${T}LT${MKL_T}MKL.txt startup_ara.sbatch $T $MKL_T true $N_RESAMPLES $N $N_SAMPLES $EPS_ANGL $EPS $N_SAMPLES_COMPOUND_EACH
sleep 0.1



