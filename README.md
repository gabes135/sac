# sac - Julia code to run the Stochastic Analytic Continuation Method
Currently only supports unconstrained sampling and the monotonic edge constrained parameterizations for fermionic spectral functions.

## Prepare G(tau) data

To run the sac_free.jl, you first need to generate a t.in file containing the covariance matrix for your G(τ) data. This is done using `make_tin.jl` in the folder `sac/process_G`, which reads in your G(τ) in the form of a full list of G(τ) bins (not bin averages). The two input files for `make_tin.jl` are `cor.dat` and `tgrid.dat` (see examples of these files in the `process_G` folder). **Note: you need to edit the inverse temperature β set in the `run` function in make_tin.jl for your specific data set**

The file `cor.dat` contains the raw G(τ) bins. The value of G(τ_i) for each bin is listed (one value per row), and each bin is seperated by "1," which is used a an indicator to seperate the bins. This will look like:
|cor.dat|
|---|
|1|
|G_1(τ_0)|
|G_1(τ_1)|
|...|
|G_1(τ_N-1)|
|G_1(β)|
|1|
|G_2(τ_0)|
|G_2(τ_1)|
|...|
|G_2(τ_N-1)|
|G_2(β)|
|1|
|...|

For proper normalization of the *finite temperature* fermionic spectral function, your G(tau) data must include both G(τ = 0) and G(τ = β).

***If you generate your G(tau) data using a zero temperature QMC method, then use make_tin_zeroT.jl, which treats the normalization of A(ω) differently. In this case, the only required τ point is G(τ = 0).***


The file `tgrid.dat` contains the τ grid that G(τ) is evaluated on. For the file, each value of τ_i is listed, one value per row. This will look like:
|tgrid.dat|
|---|
|τ_0|
|τ_1|
|...|
|τ_N-1|
|β|

Running `make_tin.jl` will generate a `t.in` file containing G(τ) averaged over all bins, the corresponding error bars, and the covariance (in imaginary time) matrix eigenvalues. You can rename this file and move into the `sac/fermion/in_files` directory.



After generating, renaming, and moving the `t.in file`, navigate to the folder corresponding to the SAC paramterization you want to us ([`sac/free`](./sac/fermion/free/README.md) or [`sac/edge`](./sac/fermion/edge/README.md)) and follow the instructions in the README.md files (follow hyperlinks above).

