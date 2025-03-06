# sac - Julia code to run the Stochastic Analytic Continuation Method
Currently only supports unconstrained sampling and the monotonic edge constrained parameterizations for fermionic spectral functions.

## Instructions for running unconstrained sampling


### Set SAC parameters 
Set the parameters for the SAC run in `in_free.in`. The inputs are:

|in_free.in|
|---|
|par: update(s) type: 1=freq, 2=freq+amp, 3=freq with non equal amps|
|N_ω, ω_0, ω_m, δω, δω_h: number of δ's, lower ω bound, upper ω bound, ω spacing, ω spacing for output spec|
|θ_0, f_anneal, f_final, a1, a2: initial sampling temp, temperature reduction factor for main anneal, temperature reduction factor for final anneal, lowest a for theta criterion, highest a for theta criterion|
|N_anneal, anneal_steps, sample_steps: Number of anneal steps, number of sweeps per anneal step in main anneal, number of sweeps per anneal step in final anneal |
|output folder: directory to write output of SAC program to|
|symm, kernel type: whether to enforce A(-omega) = A(omega) (1) are sample A(omega) across the full freq. axis (0), which type of kernel to use when converting A(ω) to G(τ) (`zeroT` or `finiteT` or `bosonic`)|

This `in_free.in` file must be in the same directory as `sac_free.jl`. See provided `in_free.in` file for example parameters.


For this verion of the program, the user must set the *kernel type*:


# Kernel Types

| Type | Transform  | Notes |
|-----------|----------------|----|
| `zeroT` | $G(\tau) = \int_{-\infty}^{\infty} e^{-\tau \omega } A(\omega)$ | $A(-\omega) = A(\omega)$ or $A(-\omega) \neq A(\omega)$  |
| `finiteT`  | $G(\tau) = \int_{-\infty}^{\infty} \frac{e^{-\tau \omega } }{1 + e^{-\beta \omega}}A(\omega)$ | $A(-\omega) = A(\omega)$ or $A(-\omega) \neq A(\omega)$  |
| `bosonic`    | $G(\tau) = \int_{-\infty}^{\infty} e^{-\tau \omega } A(\omega) =   \int_{0}^{\infty} \left(e^{-\tau \omega } + e^{-(\beta - \tau) \omega } \right) A(\omega)$ | $A(-\omega) = e^{-\beta \omega} A(\omega)$ enforced explicitly and only positive frequency axis sampled |


### Run SAC
The program is run by executing:

`julia sac_free.jl`



