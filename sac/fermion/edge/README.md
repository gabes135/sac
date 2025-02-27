# sac - Julia code to run the Stochastic Analytic Continuation Method
Currently only supports unconstrained sampling and the monotonic edge constrained parameterizations for fermionic spectral functions.

## Instructions for running the monotonic edge constrained parameterization


### Set SAC parameters 
Set the parameters for the SAC run in `in_edge.in`. The inputs are:

|in_edge.in|
|---|
|N_e, N_c: number of δ's in monotonic edge part of spectrum, number of δ's in continuum part of spectrum|
|p, A_c, A_r: power in edge singularity, (ω-ω1)^{-p}, spectral weight of continuum ∈ [0, 1), spectral weight of rightward decaying edge ∈ [0, 1]|
|ω_0, ω_m, δω, δω_h: lower ω bound, upper ω bound, ω spacing, ω spacing for output spec|
|θ_0, f_anneal, N_anneal, a: initial sampling temp, temperature reduction factor for main anneal,  Number of anneal steps, a for theta criterion for output spectrum|
|anneal_steps, sample_steps, bins:, number of sweeps per anneal step, number of sweeps in final sampling, number of bins of anneal_steps sweeps per anneal step |
|input_file, output folder: prepared filed containing QMC G(τ) data, directory to write output of SAC program to|
|fix_edge, kernel type: whether sample edge (0) or fix it to value (by setting fix_edge = #), which type of kernel to use when converting A(ω) to G(τ) ('zeroT' or 'finiteT')|
|mode: which type of edge spectrum to resolve ('single_edge', 'single_edge_out', 'double_edge_symm', 'double_edge_in'), see examples below|

This `in_edge.in` file should be in the same directory as `sac_edge.jl`. See provided `in_edge.in` file for example parameters.


For this verion of the program, the user must set both the *kernel type* (either zero for finite temperature) and the *mode* (which type of edge spectrum to resolve). The names of the modes and examples of the spectra they produce, using artificial spectra (black) and synthetic QMC data, are shown below (blue).

`single_edge` mode:

![single_edge](../../plotting/figs/t1.jpg)

`double_edge_symm` mode:

![double_edge_symm](../../plotting/figs/t2.jpg)

`double_edge_in` mode:

![double_edge_in](../../plotting/figs/t3.jpg)

### Run SAC
The program is run by executing:

`julia sac_edge.jl`


