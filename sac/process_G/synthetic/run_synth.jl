include("generate_synth.jl")
# include("process_cor.jl")

make_spec = true


name = "fermionic_asymm_peak"

kernel_type = :finiteT

folder = @sprintf "./%s" name
mkpath(folder)
println(folder)
N_b = 8000
sigma = 4

spec_type = 1#6
β = 32
τ_max = 16.
δτ = .1
grid_type = 6#1
M = 30

σ = 10.0 ^ (-sigma)
ξ = 1.0
#N_b = 100


ω0 = 1.#-1.0
ω0_n = -1.5#8.0
A0 = 0.7
ω_exp = 2.#0.75 
σ0 = 0.

A_plus = .3#1.0
A_minus = .7#0

N_G = 2#2
ωGs = [1, -.5]#, 2.5]#[1.8, 2.5]
AGs = [.2, .4]#, .5]#[.4, .368]
σGs = [.5, .5]#, 1]#[0.7, 0.5]



G_bins_file =@sprintf "%s/cor.dat" folder
spec_file = @sprintf "%s/aw.dat" folder
τ_grid_file = @sprintf "%s/tgrid.dat" folder 




synth = Synth(spec_type, β, τ_max, δτ, grid_type, M,
      σ, ξ, N_b,
      ω0, ω0_n, A0, ω_exp, σ0,
      A_plus, A_minus,
      N_G, ωGs, AGs, σGs,
      spec_file, τ_grid_file, G_bins_file,
      kernel_type)


if make_spec
      println("Generating synthetic spectrum.")

      init_weights(synth)
      println("Writing spectrum.")
      write_spec(synth)
      println("Making τ grid.")
      tau_grid(synth)
      println("Making G(τ).")
      make_G_tau(synth)
      println("Writing G(τ) bins.")
      write_Gbins(synth)
else
      println("Spectrum and bins already made.")
end
