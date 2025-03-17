# Detecting high-dimensional time-bin-entanglement in fiber-loop systems

This repository contains the data from numerical simulations presented in the paper "Detecting high-dimensional time-bin entanglement in fiber-loop systems" (https://arxiv.org/abs/2502.18336) as well as all scripts used to generate and visiualize it. The source code used in the scripts is stored seperately in the Julia package "TimeBinEncoding.jl" (https://github.com/NiklasEuler/TimeBinEncoding.jl).

Each directory in this repository contains data and scripts related to one Figure in the manuscript. The data is contained in the Data subdirectories and can be loaded and visulized by running the corresponding `*Eval.jl` PLuto-notebook files.

To run the numerical simulations again, first install "TimeBinEncoding.jl" by running `using Pkg; Pkg.add(path="https://github.com/NiklasEuler/TimeBinEncoding.jl.git")` in the Julia REPL. Each numerical simulation can then be initiated by running the local bash script `sh run_local.sh`. For completeness, the original script to run the simulation via Slurm is also included; to queue via Slurm, run `sh scheduler_ara.sh`. Last, the computation for achievable fidelities for the two SPDC pairs in contained in the Pluto notebook `MultiSPDC.jl`.
