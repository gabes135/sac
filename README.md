# sac - Julia code to run the Stochastic Analytic Continuation Method
Currently only supports unconstrained sampling and the monotonic edge constrained parameterizations for fermionic spectral functions.

## Prepare G(tau) data

To run the sac_free.jl, you first need to generate a t.in file containing the covariance matrix for your G(tau) data. This is done using make_tin.jl, which reads in a complete list of G(tau) bins. The two input files for make_tin.jl are cor.dat and tgrid.dat (see examples of these files in the exampels folder). **Note: you need to edit the inverse temperature β set in the `run` function in make_tin.jl for your specific data set**

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

For proper normalization of the fermionic spectral function, your G(tau) data must include both G(0) and G(beta).

***If you generate your G(tau) data using a zero temperature QMC method, then use make_tin_zeroT.jl, which treats the normalization of A(ω) differently***


tgrid.dat contains the tau grid that G(tau) is evaluated on. For the file, each value of tau_i is listed, one value per row. This will look like:
|tgrid.dat|
|---|
|tau_0|
|tau_1|
|...|
|tau_N-1|
|beta|

This will generate a t.in file in the `process_G` directory, which you can rename and move into the `in_files` directory.



After generating, renaming, and moving the t.in file, navigate to the folder corresponding to the SAC paramterization you want to us (`free` or `edge`) and follow the instructions in the README.md files.



