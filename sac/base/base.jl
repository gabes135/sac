module SACs

using Dates, DelimitedFiles, LinearAlgebra, Printf

export SAC, read_G!,finiteT_K, zeroT_K, bosonic_K, calc_χ2, calc_χ2_new, @swap


Base.@kwdef mutable struct SAC

    N_ω::Int64         # Number of delta functions
    ω_0::Float64       # Lower bound on deltas
    ω_m::Float64       # Upper bound on deltas
    ωi_0::Int64 = 0    # Index of lower bound on deltas
    ωi_m::Int64 = 0    # Index of upper bound on deltas
    δω::Float64        # Spacing of frequency grid for sampling 
    δω_h::Float64      # Spacing of frequency grid for output histogram of spectrum


    anneal_steps::Int64   # Number of sweeps performed in each temperature step of the annealing process
    sample_steps::Int64   # Number of sweeps performed in the final sampling stage, usually larger than anneal_steps
    θ_0::Float64          # Initial sampling temperature
    N_anneal::Int64       # Number of temperature steps in anneal
    f_anneal::Float64     # Factor θ is divided by at each annealing step: θ_{i+1} = θ_{i}/f_anneal
    tol::Float64 = 1e-3   # Used to end annealing process early, is |χ^2_{i+1} - χ^2_{i}| < tol
    a_criterion::Float64  # A value in optimal theta criterion


    G_file::String        # File containing G(τ) data and covariance matrix

    β::Float64 = 0.       # Inverse temperature
    N_τ::Int64 = 0        # Number of τ points
    norm::Float64 = 0.    # Normalization of spectral function
    τ::Array{Float64, 1} = Array{Float64}(undef, 0)   # Array of τ points
    G::Array{Float64, 1} = Array{Float64}(undef, 0)    # Average of G(τ) bins
    σ_inv::Array{Float64, 1} = Array{Float64}(undef, 0) # Inverse of covariance matrix eigenvalues
    cov::Array{Float64, 2} = Array{Float64}(undef, 0, 0) # Covariance matrix
    ω_window::Float64 = 0.  # Rough frequency scale deduced from G(τ)
    approx_ratio::Float64 = 0.  # Approximate weight of ratio of pos. and neg. peak based on ratio of G(τ=0)/G(τ=β)


    ωi_array::Array{Int64, 1} = Array{Int64}(undef, 0)   # Delta function locations stored as indices (for free and peak)
    ω_array::Array{Float64, 1} = Array{Float64}(undef, 0) # Delta function locations stored as floats (for edge)
    A_array::Array{Float64, 1} = Array{Float64}(undef, 0) # Delta function amplitudes 
    K_array::Array{Float64, 2} = Array{Float64}(undef, 0, 0) # Kernel 


    Gbar::Array{Float64, 1} = Array{Float64}(undef, 0)      # SAC G(τ) from ωi_array, A_array
    Gbar_new::Array{Float64, 1} = Array{Float64}(undef, 0)  # SAC G(τ) after update from ωi_array, A_array
    χ2::Float64 = 0.                                        # χ² value
    χ2_new::Float64 = 0.                                    # χ² value after update
    χ2_min::Float64 = 0.                                    # Minimum χ² during anneal

    χ2_anneal::Array{Float64, 1} = Array{Float64}(undef, 0)  # χ² vs. Theta
    sampled_χ2::Float64 = 0.                                # Accumulated χ² during sampling


    anneal_file::String   # Output file for annealing results
    sample_file::String   # Output file for info on final sampling
    output_folder::String # Folder to write these files to

    kernel_type::Symbol    # finiteT or zeroT or bosonic
            
    symm::Int64           # For fermionic spectral functions only!
                          # Whether or not to impose particle-hole symmetry for Hubbard-like models
                          # If set to true, A(-ω) = A(ω) and only the positive part of the spectrum is sampled


end


function read_G!(self)
    G_file = self.G_file
    t_in = readdlm(G_file)
    
    self.β, self.N_τ, n_bins, self.norm = t_in[1,:]
    
    N_τ = self.N_τ
    
    cov = zeros((N_τ, N_τ))
    
    τ, G, σ = t_in[2:N_τ+1,1], t_in[2:N_τ+1,2], t_in[2:N_τ+1,4]
    σ_inv = 1. ./ σ
    
    cov_start = N_τ+2
    
    for j=0:N_τ-1
        @assert j+1 == t_in[cov_start + j*(N_τ+1),1]
        cov[:, j+1] = t_in[cov_start + j*(N_τ+1) + 1:cov_start + j*(N_τ+1) + N_τ,1]
    end
    
    G_half = G[τ .<= self.β ÷ 4]
    τ_half = τ[τ .<= self.β ÷ 4]
    
    ω_window = log(1/G_half[end])/τ_half[end]


    if self.symm == 1 || self.kernel_type == :bosonic || maximum(τ) <= self.β ÷ 2
        # neg. axis not sampled so all peak weight on pos. axis
        approx_ratio = Inf
    else
        approx_ratio = G[1]/G[N_τ]
    end
    
    G_D = transpose(cov) * G
    
    self.τ = τ
    self.G = G_D
    self.σ_inv = σ_inv
    self.cov = cov
    self.ω_window = ω_window
    self.approx_ratio = approx_ratio
    
end


macro swap(x,y)
   quote
      local tmp = $(esc(x))
      $(esc(x)) = $(esc(y))
      $(esc(y)) = tmp
    end
end



function finiteT_K(ω, τ, β)
    num = -ω * τ
    denom = -β * ω

    max_term = max(num, denom, 0) 
 
    num_stab = exp(num - max_term)
    denom_stab = exp(-max_term) + exp(denom - max_term)

    return num_stab / denom_stab
end


function zeroT_K(ω, τ, β)
    return exp(-ω * τ)
end

function bosonic_K(ω, τ, β)
    return (exp(-ω * τ) + exp(-ω * (β - τ)))/(1 + exp(-β * ω))
end

# χ^2 = Σ (G_bar - G) / σ^2
function calc_χ2(self)
    χ2 = 0
    @inbounds @simd for i=1:self.N_τ
        χ2 += ((self.Gbar[i] - self.G[i]) * self.σ_inv[i]) ^ 2
    end
    return χ2
end

function calc_χ2_new(self)
    χ2 = 0
    @inbounds @simd for i=1:self.N_τ
        χ2 += ((self.Gbar_new[i] - self.G[i]) * self.σ_inv[i]) ^ 2
    end
    return χ2
end


end 