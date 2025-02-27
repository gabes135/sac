include("generate_synth.jl")
# include("process_cor.jl")

make_spec = true
make_tin = false

name = "double_edge_out"

folder = @sprintf "./%s" name
mkpath(folder)
println(folder)
N_b = 8000
sigma = 3

spec_type = 3#6
β = 32#128
τ_max = 16.#8.
δτ = .1
grid_type = 6#1
M = 30

σ = 10.0 ^ (-sigma)
ξ = 1.0
#N_b = 100


ω0 = 1.#-1.0
ω0_n = -1.5#8.0
A0 = 0.5
ω_exp = 2.#0.75 
σ0 = .25

A_plus = 0.7#1.0
A_minus = 0.3#0

N_G = 1#2
ωGs = [0.]#, 2.5]#[1.8, 2.5]
AGs = [0.]#, .5]#[.4, .368]
σGs = [0.]#, 1]#[0.7, 0.5]



G_bins_file =@sprintf "%s/cor.dat" folder
spec_file = @sprintf "%s/aw.dat" folder
τ_grid_file = @sprintf "%s/tgrid.dat" folder 




synth = Synth(spec_type, β, τ_max, δτ, grid_type, M,
      σ, ξ, N_b,
      ω0, ω0_n, A0, ω_exp, σ0,
      A_plus, A_minus,
      N_G, ωGs, AGs, σGs,
      spec_file, τ_grid_file, G_bins_file)


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
# println()
# if make_tin
#       τ_min = 0.0
#       r_b = 1
#       τ_freq = 1
#       N_boot = 10000

#       cor_file = G_bins_file
#       log_file = ""
#       G_norm_file = ""
#       G_sigma_file = @sprintf "%s/G.dat" folder
#       out_file = @sprintf "%s/t.in" folder

#       println("Processing synthetic spectrum.")

#       tcor = TimeCor(β, τ_min, r_b, τ_freq, N_boot,
#                       0, N_b, 
#                       τ_grid_file, cor_file, log_file,
#                       G_norm_file, G_sigma_file, out_file)

#       println("Reading G(τ) bins.")
#       read_data(tcor)
#       println("Computing G(τ) average and normalization.")
#       compute_means(tcor)
#       println("Constructing convariance matrix.")
#       cov_matrix(tcor)
# end
