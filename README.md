# sac - Julia code to run the Stochastic Analytic Continuation Method
Currently only supports unconstrained sampling for fermionic spectral functions.

## Step 1: Prepare G(tau) data

To run the sac_free.jl, you first need to generate a t.in file containing the covariance matrix for your G(tau) data. This is done using make_tin.jl, which reads in a complete list of G(tau) bins. The two input files for make_tin.jl are cor.dat and tgrid.dat (see examples of these files in the exampels folder). **Note: you need to edit the inverse temperature set in the `run` function in make_tin.jl for your specific data set**

cor.dat contains the G(tau) bins. The value of G(tau_i) for each bin is listed (one value per row), and each bin is seperated by "1," which is just usd a an indicator to seperate the bins. This will look like:
|cor.dat|
|---|
|1|
|G_1(tau_0)|
|G_1(tau_1)|
|...|
|G_1(tau_N-1)|
|G_1(beta)|
|1|
|G_2(tau_0)|
|G_2(tau_1)|
|...|
|G_2(tau_N-1)|
|G_2(beta)|
|1|
|...|

For proper normalization of the Fermionic spectral function, your G(tau) data must include both G(0) and G(beta).


tgrid.dat contains the tau grid that G(tau) is evaluated on. For the file, each value of tau_i is listed, one value per row. This will look like:
|tgrid.dat|
|---|
|tau_0|
|tau_1|
|...|
|tau_N-1|
|beta|

This t.in file should be in the same directory as sac_free.jl.


## Step 2: Set SAC parameters 
After generating the t.in file, set the parameters for the SAC run in in_free.in. The inputs are (see provided example):

|in_free.in|
|---|
|par: update(s) type: 1=freq, 2=freq+amp, 3=freq with non equal amps|
|N_ω, ω_0, ω_m, δω, δω_h: number of δ's, lower ω bound, upper ω bound, ω spacing, ω spacing for output spec|
|θ_0, f_anneal, f_final, a1, a2: initial sampling temp, temperature reduction factor for main anneal, temperature reduction factor for final anneal, lowest a for theta criterion, highest a for theta criterion|
|N_anneal, anneal_steps, sample_steps: Number of anneal steps, number of sweeps per anneal step in main anneal, number of sweeps per anneal step in final anneal |
|output folder: directory to write output of SAC program to|
|phs: whether to enforce A(-omega) = A(omega) (1) are let A(-omega) differ from A(omega) (0)|

This in_free.in file should be in the same directory as sac_free.jl


## Step 3: Run SAC
The program is run by executing:

`julia sac_free.jl`



