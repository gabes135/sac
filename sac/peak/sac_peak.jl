#FERMIONIC

using DelimitedFiles, LinearAlgebra, Statistics, Printf, Dates

module SACs
Base.@kwdef mutable struct SAC
    
    N_ω::Int64 # Number of delta functions in continuum
    N_p::Int64 # Number of delta functions in leadging peak
    ω_0::Float64 = 0. # Lower bound on δ's
    ω_m::Float64 # Upper bound on δ's
    ωi_0::Int64 = 0 # Index of lower bound on δ's
    ωi_m::Int64 = 0 # Index of upper bound on δ's
    ωi_pp::Float64 = 0 # Index of right-most (max) peak delta on positive freq. axis
    ωi_pc::Int64 = 0 # Index of left-most (min) continuum delta on positive freq. axis
    ωi_np::Float64 = 0 # Same as above but for negative freq. axis
    ωi_nc::Int64 = 0 # Same as above but for negative freq. axis
    δω::Float64 # Spacing of frequency grid for sampling 
    δω_h::Float64 # Spacing of frequency grid for output histogram of spectrum
    A_0::Float64 # Weight of leading peak feature(s) (combined weight on pos. and neg. freq. axis)

    peak_p::UnitRange{Int64} = 0:0 # Indices of pos. peak δ's in ωi_array
    cont_p::UnitRange{Int64} = 0:0 # Indices of pos. cont. δ's in ωi_array
    peak_n::UnitRange{Int64} = 0:0 # Indices of neg. peak δ's in ωi_array
    cont_n::UnitRange{Int64} = 0:0 # Indices of neg. cont. δ's in ωi_array

    anneal_steps::Int64 # Number of sweeps preformed in each temperture step of the annealing processs
    sample_steps::Int64 # Number of sweeps preformed in the final sampling stage, usually larger than anneal_steps
    θ_0::Float64 # Initial sampling temperature
    N_anneal::Int64 # Number of temperature steps in anneal
    f_anneal::Float64 # Factor θ is divided by at each annealing steps: θ_{i+1} = θ_{i}/f_anneal
    
    tol::Float64 = 1e-3 # Used to end annealing process early, is |χ^2_{i+1} - χ^2_{i}| < tol
    
    a_criterion::Float64 # a value in optimal theta criterion
  
    G_file::String # File containing G(tau) data and covariance matrix

    β::Float64 = 0. # Inverse temperature
    N_τ::Int64 = 0 # Number of τ points
    norm::Float64 = 0. # Normalization of spectral function
    τ::Array{Float64, 1} = Array{Float64}(undef, 0) # Array of τ points
    G::Array{Float64, 1} = Array{Float64}(undef, 0) # Average of G(τ) bins
    σ_inv::Array{Float64, 1} = Array{Float64}(undef, 0) # Inverse of covariance matrix eigenvalues
    cov::Array{Float64, 2} = Array{Float64}(undef, 0, 0) # Covariance matrix
    ω_window::Float64 = 0. # Rough frequency scale deduced from G(τ)
    approx_ratio::Float64 = 0. # Approximate weight of ratio of pos. and neg. peak based on ratio of G(τ=0)/G(τ=β)
    aprox_edges::Array{Float64, 1} = zeros(Float64, 2) 

    ωi_array::Array{Int64, 1} = Array{Int64}(undef, 0)  # Delta function locations stored as integers
    A_array::Array{Float64, 1} = Array{Int64}(undef, 0) # Delta function amplitudes 
    Kp_array::Array{Float64, 2} = Array{Int64}(undef, 0, 0) # Kernel for pos. freq. contribution
    Kn_array::Array{Float64, 2} = Array{Int64}(undef, 0, 0) # Kernel for neg. freq. contribution

    Gbar::Array{Float64, 1} = Array{Float64}(undef, 0) #SAC G(τ) from ωi_array, A_array
    Gbar_new::Array{Float64, 1} = Array{Float64}(undef, 0) #SAC G(τ) after updatefrom ωi_array, A_array
    χ2::Float64 = 0. 
    χ2_new::Float64 = 0. # χ^2 after update
    χ2_min::Float64 = 0. # Minimum χ^2 during anneal
    χ2_anneal::Array{Float64, 1} = Array{Float64}(undef, 0) # χ^2 vs. Theta

    accept_rates::Array{Float64, 1} =  zeros(Float64, 11) # Acceptance rates for various updates
    update_windows::Array{Float64, 1} =  zeros(Float64, 11) # Frequency and ampltiude windows used for updates

    sampled_pspec::Array{Float64, 2} = Array{Float64}(undef, 0, 0) # Accumlated peak spectrum during sampling
    sampled_cspec::Array{Float64, 2} = Array{Float64}(undef, 0, 0) # Accumlated cont. spectrum during sampling
    sampled_χ2::Float64 = 0. # Accumulated χ^2 during sampling

    anneal_file::String # Output file for annealing results
    sample_file::String # Output file for info on final sampling
    accept_rate_file::String # Output file for acceptance rates
    output_folder::String # Folder to write these files to

    edge_res::Array{Float64, 2} = zeros(3, 2) # Sampled location(s) of peak edge(s)

    kernel_type::Symbol # finiteT or zeroT or bosonic

    symm::Int64 # For fermionic spectral functions only!
                # Whether or not to impose particle-hole symettry for Hubbard like models
                # If set to true, A(-ω) = A(ω) and only the positive part of the spectrum is sampled
     
    indiv_update::Bool = false
end
end
SAC = SACs.SAC

###########################################################################################
## Intialization Functions
###########################################################################################


# Read the input file "G_file"
# This file should be prepared by "make_tin.jl"
function read_G!(self::SAC)
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
    # approx_ratio = 1.
    G_D = transpose(cov) * G
    
    self.τ = τ
    self.G = G_D
    self.σ_inv = σ_inv
    self.cov = cov
    self.ω_window = ω_window
    self.approx_ratio = approx_ratio
    
end

# Function to determine position of peak position for pos. (pn=1) or neg. (pn=2) axis
function peak_bound(self::SAC, pn::Int64)
    if pn == 1
        return maximum(Int64, self.ωi_array[self.peak_p])
    else
        return maximum(Int64, self.ωi_array[self.peak_n])
    end
end

function initialize!(self::SAC)
    
    ωi_0 = floor(Int64, self.ω_0/self.δω) # Converting ω_0 to integer
    ωi_m = ceil(Int64, self.ω_m/self.δω)
    
    ω_window = self.ω_window/self.δω
    
    # Array containing δ locations and amplitudes 
    # Order is positive peak, pos. cont., neg. peak, neg cont.
    ωi_array = zeros(Int64, 2*(self.N_ω + self.N_p)) 
    A_array = zeros(Float64, 2*(self.N_ω + self.N_p))

    # Accumated spectra for peak and cont., second index is pos./neg. axis
    self.sampled_pspec = zeros(Float64, (ωi_m+1), 2)
    self.sampled_cspec = zeros(Float64, (ωi_m+1), 2)
    
    # Indices of ωi_array and A_array corresponding to the four different contributions to A(ω)
    self.peak_p = 1:self.N_p
    self.cont_p = (1+self.N_p):(self.N_ω+self.N_p)
    self.peak_n = (1+self.N_ω+self.N_p):(self.N_ω+2*self.N_p)
    self.cont_n = (1+self.N_ω+2*self.N_p):(2*self.N_ω+2*self.N_p)

    # Initial guess for location and amplitude of pos./neg. peaks
    A0_p = self.A_0/(1 + (1/self.approx_ratio))
    A0_n = self.A_0/(1 + self.approx_ratio)



    ωi_array[self.peak_p] .= floor(Int64, ω_window)#floor(Int64, 1/self.δω)#
    A_array[self.peak_p] .= A0_p/self.N_p # distributed evenly amongst the peak δ's

    ωi_array[self.peak_n] .= floor(Int64, ω_window)#floor(Int64, .5/self.δω)#floor(Int64, ω_window)
    A_array[self.peak_n] .= A0_n/self.N_p # distributed evenly amongst the peak δ's
  
    # Right-most edge of peak feature
    ωi_pp = maximum(ωi_array[self.peak_p])
    ωi_np = maximum(ωi_array[self.peak_n])

    # First set initial loc of cont. to just above edge
    ωi_array[self.cont_p] .= floor(Int64, ωi_pp*1.0)
    ωi_array[self.cont_n] .= floor(Int64, ωi_np*1.0)

    # Now spread them out with spacing proportionall to ω_window
    ωi_array[self.cont_p] .+= floor.(Int64, (ω_window/self.N_ω) .* (1:self.N_ω))
    ωi_array[self.cont_n] .+= floor.(Int64, (ω_window/self.N_ω) .* (1:self.N_ω))
    
    # Left-most edge of cont.
    ωi_pc = minimum(ωi_array[self.cont_p]) 
    ωi_nc = minimum(ωi_array[self.cont_n]) 

    # Distribute remaining weight evenly amongst all continu (equal ration between pos./neg.)
    Ac = (1-self.A_0)/2
    Ac_p = (1-self.A_0)/(1 + (1/self.approx_ratio))
    Ac_n = (1-self.A_0)/(1 + self.approx_ratio)

    A_array[self.cont_p] .= (1:self.N_ω)
    A_array[self.cont_n] .= (1:self.N_ω)
    
    A_array[self.cont_p] .*= Ac_p/sum(1:self.N_ω) 
    A_array[self.cont_n] .*= Ac_n/sum(1:self.N_ω)

    # If symm, then only sample positive freq. axis
    if self.symm == 1 || self.kernel_type == :bosonic
        A_array[self.peak_p] .= self.A_0/self.N_p
        A_array[self.cont_p] .= (1-self.A_0)/self.N_ω
        A_array[self.peak_n] .= 0
        A_array[self.cont_n] .= 0
        if self.symm == 1
            A_array ./= (2*sum(A_array)) # If symm, then only positive side of frequency axis is sampled,
                                         # and is copied to negative frequency axis.
                                         # Need to divide amplitude by two to maintain normaliztion
                                         # Σ A_i = 1
        end

        ωi_array[self.peak_n] .= 0
        ωi_array[self.cont_n] .= 0
    else
        A_array ./= sum(A_array)
    end
    
    # println([A0_p, A0_n])
    # println([ωi_pp, ωi_np].* self.δω)
    # println([Ac_p, Ac_n])
    # Initializing Kernel


    if self.kernel_type == :finiteT
        K_function = finiteT_K
    elseif self.kernel_type == :zeroT
        K_function = zeroT_K
    elseif self.kernel_type == :bosonic
        K_function = bosonic_K
    else
        throw("Invalid Kernel type.")
    end
    
    Kp_array = zeros(Float64, (self.N_τ, ωi_m+1)) 
    Kn_array = zeros(Float64, (self.N_τ, ωi_m+1)) 
    K_vals = zeros(Float64, self.N_τ)
   
    for i=0:ωi_m
        ω = (i+.5) * self.δω
        K_vals .= K_function.(ω, self.τ, self.β)

        # If symm, then kernel is adjusted to account for S(-ω) = S(ω) 
        # and Kn_array is left empty 
        if self.symm == 1  
            K_vals .+=  K_function.(-ω, self.τ, self.β)

            Kp_array[:, i+1] .= K_vals
            Kp_array[:, i+1] .= transpose(self.cov) * Kp_array[:, i+1] # Transform into eigenbasis of covariance matrix

        else
            Kp_array[:, i+1] .= K_vals
            Kp_array[:, i+1] .= transpose(self.cov) * Kp_array[:, i+1] # Transform into eigenbasis of covariance matrix

            # If bosonic, relation between pos./neg. spectrum imposed in kernel and only pos. axis sampled
            if self.kernel_type != :bosonic
                K_vals .= K_function.(-ω, self.τ, self.β)
                Kn_array[:, i+1] .= K_vals
                Kn_array[:, i+1] .= transpose(self.cov) * Kn_array[:, i+1]
            end

        end
    end

    # println(ωi_array)

    self.ωi_array, self.A_array, self.Kp_array, self.Kn_array = ωi_array, A_array, Kp_array, Kn_array
    self.ω_window, self.ωi_m = ω_window, ωi_m
    self.ωi_pp, self.ωi_pc, self.ωi_np, self.ωi_nc = ωi_pp, ωi_pc, ωi_np, ωi_nc
    self.χ2_anneal = zeros(self.N_anneal)
    
end

# Funtion to calculate kernel that can handle large negative ω
# K(ω, τ) = exp(-τω) / (1 + exp(-βω))
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


function init_outfiles(self)

    
    # Last four columns are peak weight on pos. axis, peak weight on neg. axis, cont. on pos., and cont. on neg.
    f = open(self.anneal_file, "w")
    writedlm(f, [["i", "theta", "chi2_min", "chi2_avg", "edge_p", "edge_n", "Ap_p", "Ap_n", "Ac_p", "Ac_n"]], ',')
    close(f)
    
    # Same as above but now also write fluctuations in edge locs
    f = open(self.sample_file, "w")
    writedlm(f, [["i", "a", "theta", "chi2_min", "chi2_avg", "edge_p", "edge_n", "sigma_edge_p", "sigma_edge_n",
                    "Ap_p", "Ap_n", "Ac_p", "Ac_n"]], ',')
    close(f)
   
   
    #acceptance rates:
        # cont. freq. (single, double, triple) for pos. and neg., and their update windows
        # peak freq. for pos. and neg., and their update windows
        # amplitude transfer for cont. and peak., and their update windows
    f = open(self.accept_rate_file, "w")
    writedlm(f, [["i", "ar_pcont1", "ar_pcont2", "ar_pcont3",
                       "ar_ncont1", "ar_ncont2", "ar_ncont3",
                       "∆w_pcont1", "∆w_pcont2", "∆w_ncont1", "∆w_ncont2",
                       "ar_ppeak", "ar_npeak",
                       "∆w_ppeak", "∆w_npeak",
                       "ar_wAc", "ar_wAp", "ar_wApc",
                       "∆w_Ac", "∆w_Ap", "∆w_Apc"]], ',')
    close(f)

 
    f = open(string(self.output_folder, "/log.txt"), "w")
    close(f)
    
end


# Calculate G(τ) for delta config, adds contributions from all four features
function calc_Gbar(self::SAC)
    Gbar = (self.Kp_array[:, self.ωi_array[self.peak_p] .+ 1] * self.A_array[self.peak_p] .+ 
            self.Kp_array[:, self.ωi_array[self.cont_p] .+ 1] * self.A_array[self.cont_p] .+ 
            self.Kn_array[:, self.ωi_array[self.peak_n] .+ 1] * self.A_array[self.peak_n] .+ 
            self.Kn_array[:, self.ωi_array[self.cont_n] .+ 1] * self.A_array[self.cont_n])
    return Gbar 
end


# Calculate χ2 for current config
function calc_χ2(self::SAC)
    χ2 = 0
    @inbounds @simd for i=1:self.N_τ
        χ2 += ((self.Gbar[i] - self.G[i]) * self.σ_inv[i]) ^ 2
    end
    return χ2
end

# Calculate χ2 for updated config, corresponding to Gbar_new
function calc_χ2_new(self::SAC)
    χ2 = 0
    @inbounds @simd for i=1:self.N_τ
        χ2 += ((self.Gbar_new[i] - self.G[i]) * self.σ_inv[i]) ^ 2
    end
    return χ2
end


###########################################################################################
## Update Functions
###########################################################################################
# pn argument denotes which freq. axis to sample (pn=1 - pos., pn = 2 - neg.)


# Update freq. of one δ in cont.
function single_ω_move(self::SAC, θ::Float64, pn::Int64)
    
    accept_rate = 0
    
    update_window = self.update_windows[1 + (pn-1)*7] # window 1 or 8

    if pn == 1
        ind = self.cont_p
        ω_b = self.ωi_pp # cont. deltas can't go post peak edge
        K_array = self.Kp_array
    else
        ind = self.cont_n
        ω_b = self.ωi_np
        K_array = self.Kn_array
    end

    ω_1 = -1
    for i=1:self.N_ω
        ωi = rand(ind)
        δωi = 1 + Int64(floor(rand()*update_window))
        
        if rand() < 0.5
            ω_1 = self.ωi_array[ωi] + δωi
        else
            ω_1 = self.ωi_array[ωi] - δωi
        end
        
        if ω_1 <= floor(Int64, ω_b)
            continue
        elseif ω_1 > self.ωi_m
            continue
        end
        
        A = self.A_array[ωi]
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
              (A * (K_array[j, ω_1 + 1] - K_array[j, self.ωi_array[ωi] + 1])))
        end
        
        χ2_new = calc_χ2_new(self)
        
        @fastmath P = exp((self.χ2 - χ2_new)/(2*θ))
     
        if rand() <= P
            self.ωi_array[ωi] = ω_1
            
            @views self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new
            
            if χ2_new < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1
        end
  
    end

    # Update continum left edge after updating
    if pn == 1
        self.ωi_pc = minimum(self.ωi_array[self.cont_p])
    else
        self.ωi_nc = minimum(self.ωi_array[self.cont_n])
    end

    return accept_rate/self.N_ω
end

# Update freq. of two δ's in cont.
function double_ω_move(self::SAC, θ::Float64,  pn::Int64)

    accept_rate = 0
    update_window = self.update_windows[2 + (pn-1)*7] # window 2 or 9

    if pn == 1
        ind = self.cont_p
        ω_b = self.ωi_pp
        K_array = self.Kp_array
    else
        ind = self.cont_n
        ω_b = self.ωi_np
        K_array = self.Kn_array
    end

    ω_1 = -1
    ω_2 = -1

    for i=1:(self.N_ω ÷ 2)
        ωi_1 = rand(ind)
        ωi_2 = ωi_1
        while ωi_2 == ωi_1
            ωi_2 = rand(ind)
        end
        
        δωi = 1 + Int64(floor(rand()*update_window))
        
        if rand() < 0.5
            ω_1 = self.ωi_array[ωi_1] + δωi
            ω_2 = self.ωi_array[ωi_2] - δωi
        else
            ω_1 = self.ωi_array[ωi_1] - δωi
            ω_2 = self.ωi_array[ωi_2] + δωi
        end
        
        if ω_1 <= floor(Int64, ω_b)
            continue
        elseif ω_2 <= floor(Int64, ω_b)
            continue
        elseif ω_1 > self.ωi_m
            continue
        elseif ω_2 > self.ωi_m
            continue
        end
        
        A1 = self.A_array[ωi_1]
        A2 = self.A_array[ωi_2]
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
                A1 * (K_array[j, ω_1 + 1] - K_array[j, self.ωi_array[ωi_1] + 1]) +
                A2 * (K_array[j, ω_2 + 1] - K_array[j, self.ωi_array[ωi_2] + 1]))
        end

        
        χ2_new  = calc_χ2_new(self)
        
        @fastmath P = exp((self.χ2 - χ2_new)/(2*θ))
        
        if rand() <= P
            self.ωi_array[ωi_1] = ω_1
            self.ωi_array[ωi_2] = ω_2
            
            @views self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new
            
            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 2
        end
        
    end

    if pn == 1
        self.ωi_pc = minimum(self.ωi_array[self.cont_p])
    else
        self.ωi_nc = minimum(self.ωi_array[self.cont_n])
    end

    return accept_rate/(self.N_ω ÷ 2)
    
    
end

# Update freq. of three δ's in cont.
function triple_ω_move(self::SAC, θ::Float64, pn::Int64)

    accept_rate = 0

    if pn == 1
        ind = self.cont_p
        ω_b = self.ωi_pp
        K_array = self.Kp_array
    else
        ind = self.cont_n
        ω_b = self.ωi_np
        K_array = self.Kn_array
    end


    ω_1 = -1 
    ω_2 = -1
    ω_3 = -1

    
    for i=1:(self.N_ω ÷ 3)
        ωi_1 = rand(ind)
        ωi_2 = ωi_1
        while ωi_2 == ωi_1
            ωi_2 =  rand(ind)
        end
        ωi_3 = ωi_2
        while ωi_3 == ωi_1 || ωi_3 == ωi_2
            ωi_3 =  rand(ind)
        end
        
        δωi = (self.ωi_array[ωi_2] + self.ωi_array[ωi_3] - (2 * self.ωi_array[ωi_1])) ÷ 3
        
        ω_1 = self.ωi_array[ωi_1] + 2*δωi
        ω_2 = self.ωi_array[ωi_2] - δωi
        ω_3 = self.ωi_array[ωi_3] - δωi
        
        if ω_1 <= floor(Int64, ω_b)
            continue
        elseif ω_2 <= floor(Int64, ω_b)
            continue
        elseif ω_3 <= floor(Int64, ω_b)
            continue
        elseif ω_1 > self.ωi_m
            continue
        elseif ω_2 > self.ωi_m
            continue
        elseif ω_3 > self.ωi_m
            continue
        end
        
        A1 = self.A_array[ωi_1]
        A2 = self.A_array[ωi_2]
        A3 = self.A_array[ωi_3]
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
                A1 * (K_array[j, ω_1 + 1] - K_array[j, self.ωi_array[ωi_1] + 1]) +
                A2 * (K_array[j, ω_2 + 1] - K_array[j, self.ωi_array[ωi_2] + 1]) +
                A3 * (K_array[j, ω_3 + 1] - K_array[j, self.ωi_array[ωi_3] + 1]))
        end

        χ2_new  = calc_χ2_new(self)
        
        @fastmath P = exp((self.χ2 - χ2_new)/(2*θ))
        
        if rand() <= P
            self.ωi_array[ωi_1] = ω_1
            self.ωi_array[ωi_2] = ω_2
            self.ωi_array[ωi_3] = ω_3
            
            @views self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new
            
            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1
        end
        
    end

    if pn == 1
        self.ωi_pc = minimum(self.ωi_array[self.cont_p])
    else
        self.ωi_nc = minimum(self.ωi_array[self.cont_n])
    end

    return accept_rate/(self.N_ω ÷ 3) 
    
  
end

# Update freq. of one δ's in peak
function single_ω_pmove(self::SAC, θ::Float64, pn::Int64)
    
    accept_rate = 0

    

     if pn == 1
        ind = self.peak_p
        ω_b = self.ωi_pc # peak δ's can't go past continuum edge
        K_array = self.Kp_array
        update_window = self.update_windows[4]
    else
        ind = self.peak_n
        ω_b = self.ωi_nc
        K_array = self.Kn_array
        update_window = self.update_windows[5]
    end

    ω_1 = -1

    for i=1:10*self.N_p
        ωi = rand(ind)
        δωi = 1 + Int64(floor(rand()*update_window))
        
        if rand() < 0.5
            ω_1 = self.ωi_array[ωi] + δωi 
        else
            ω_1 = self.ωi_array[ωi] - δωi
        end

       
        if ω_1 < self.ωi_0
            continue
        elseif ω_1 > ω_b
            continue
        end


        
        A = self.A_array[ωi]
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
              (A * (K_array[j, ω_1 + 1] - K_array[j, self.ωi_array[ωi] + 1])))
        end

        χ2_new = calc_χ2_new(self)
        
        @fastmath P = exp((self.χ2 - χ2_new)/(2*θ))
  
        if rand() <= P
            self.ωi_array[ωi] = ω_1
           
            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new
            
            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1
        end
  
    end

 
    # Update peak bound after updates
    if pn == 1
        self.ωi_pp = peak_bound(self, 1)
    else
        self.ωi_np = peak_bound(self, 2)
    end
            
    return accept_rate/(10*self.N_p)
end



# Update that transfers weight between pos. and neg. cont.'s while preserving first moment
# of two updated δ's'
function cont_Aω_transfer(self::SAC, θ::Float64)

    accept_rate = 0.
    update_window = self.update_windows[6]


    # Ai_p, Ai_p, ω_p, ω_n, ω_p_prime, ω_n_prime = 0, 0, 0, 0, 0, 0
    # A_p, A_n, A_p_prime, A_n_prime = 0., 0., 0., 0.
    # δ_A, δωi = 0., 0.

   
    for i=1:(self.N_ω ÷ 2)
        Ai_p = rand(self.cont_p)
        Ai_n = rand(self.cont_n)

        ω_p = self.ωi_array[Ai_p]
        ω_n = self.ωi_array[Ai_n]
         
        δωi = 1 + floor(Int64, rand()*update_window)

        if rand() < 0.5
            δωi *= -1
        end
        
        ω_p_prime = ω_p + δωi
        ω_n_prime = ω_n - δωi
      

        if ω_p_prime > self.ωi_m || ω_p_prime < self.ωi_0 || ω_p_prime < self.ωi_pp
            continue
        end

        if ω_n_prime > self.ωi_m || ω_n_prime < self.ωi_0 || ω_n_prime < self.ωi_pp
            continue
        end

        A_p = self.A_array[Ai_p] 
        A_n = self.A_array[Ai_n] 

        δ_A = (δωi * (A_n - A_p) / (2*δωi + self.ωi_array[Ai_p] - (-self.ωi_array[Ai_n])))  #To conserve first moment


        A_p_prime = A_p + δ_A
        A_n_prime = A_n - δ_A

        if A_p_prime < 0
            continue
        elseif A_n_prime < 0
            continue
        end

        
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j]
                    + (A_p_prime * self.Kp_array[j, ω_p_prime + 1])
                    + (A_n_prime * self.Kn_array[j, ω_n_prime + 1])
                    - (A_p * self.Kp_array[j, ω_p + 1]) 
                    - (A_n * self.Kn_array[j, ω_n + 1]))
        end

        χ2_new = calc_χ2_new(self)
        
        P = exp((self.χ2 - χ2_new)/(2*θ))


        if rand() <= P
            self.ωi_array[Ai_p] = ω_p_prime
            self.ωi_array[Ai_n] = ω_n_prime
            self.A_array[Ai_p] = A_p_prime
            self.A_array[Ai_n] = A_n_prime

            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new

            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1.
        end
    end

    self.ωi_pc = minimum(self.ωi_array[self.cont_p])
    self.ωi_nc = minimum(self.ωi_array[self.cont_n])
    
    return accept_rate / (self.N_ω ÷ 2)
    
end


# Update that transfers weight between pos. and neg. peaks's while preserving first moment
# of two updated δ's
function peak_Aω_transfer(self::SAC, θ::Float64)

    accept_rate = 0.
    update_window = self.update_windows[7]

    for i=1:(10) 
        Ai_p = rand(self.peak_p)
        Ai_n = rand(self.peak_n)

        ω_p = self.ωi_array[Ai_p]
        ω_n = self.ωi_array[Ai_n]
         
        δωi = 1 + floor(Int64, rand()*update_window)


        if rand() < 0.5
            δωi *= -1
        end

        ω_p_prime = ω_p + δωi
        ω_n_prime = ω_n - δωi

        if ω_p_prime < self.ωi_0 || ω_p_prime > self.ωi_pc || ω_p_prime > self.ωi_m
            continue
        end
        if ω_n_prime < self.ωi_0 || ω_n_prime > self.ωi_nc || ω_n_prime > self.ωi_m
            continue
        end

        A_p = self.A_array[Ai_p] 
        A_n = self.A_array[Ai_n]


        δ_A = (δωi * (A_n - A_p) / (2*δωi + self.ωi_array[Ai_p] - (-self.ωi_array[Ai_n])))  #To conserve first moment

     
        A_p_prime = A_p + δ_A
        A_n_prime = A_n - δ_A

        if A_p_prime < 0
            continue
            # return accept_rate
        elseif A_n_prime < 0
            continue
            # return accept_rate
        end

        
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j]
                    + (A_p_prime * self.Kp_array[j, ω_p_prime + 1])
                    + (A_n_prime * self.Kn_array[j, ω_n_prime + 1])
                    - (A_p * self.Kp_array[j, ω_p + 1]) 
                    - (A_n * self.Kn_array[j, ω_n + 1]))
        end

        χ2_new = calc_χ2_new(self)
        
        P = exp((self.χ2 - χ2_new)/(2*θ))


        if rand() <= P
            self.ωi_array[Ai_p] = ω_p_prime
            self.ωi_array[Ai_n] = ω_n_prime
            self.A_array[Ai_p] = A_p_prime
            self.A_array[Ai_n] = A_n_prime

            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new

            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1.
        end
    end


    self.ωi_pp = peak_bound(self, 1)    
    self.ωi_np = peak_bound(self, 2)

    return accept_rate/10
end


# Update that transfers weight between pos. and neg. peaks and continua at same time, while preserving first moment
# of four updated δ's
function Aω_transfer(self::SAC, θ::Float64)

    accept_rate = 0.
    c_update_window = self.update_windows[11]
    p_update_window = self.update_windows[11]


   
    for i=1:(self.N_ω ÷ 2)
        Ai_pc = rand(self.cont_p)
        Ai_nc = rand(self.cont_n)

        ω_pc = self.ωi_array[Ai_pc]
        ω_nc = self.ωi_array[Ai_nc]

        Ai_pp = rand(self.peak_p)
        Ai_np = rand(self.peak_n)

        ω_pp = self.ωi_array[Ai_pp]
        ω_np = self.ωi_array[Ai_np]
         
        δωi_c = 1 + floor(Int64, rand()*c_update_window)
        δωi_p = 1 + floor(Int64, rand()*p_update_window)

        if rand() < 0.5
            δωi_c *= -1
        end

        if rand() < 0.5
            δωi_p *= -1
        end
        
       
        ω_pc_prime = ω_pc + δωi_c
        ω_nc_prime = ω_nc - δωi_c
        
        if ω_pc_prime > self.ωi_m || ω_pc_prime < self.ωi_pp || ω_pc_prime < self.ωi_0
            continue
        end
        if ω_nc_prime > self.ωi_m || ω_nc_prime < self.ωi_np || ω_nc_prime < self.ωi_0
            continue
        end

        ω_pp_prime = ω_pp + δωi_p
        ω_np_prime = ω_np - δωi_p
        

        if ω_pp_prime > self.ωi_pc || ω_pp_prime > self.ωi_m || ω_pp_prime < self.ωi_0
            continue
        end
        if ω_np_prime > self.ωi_nc || ω_np_prime > self.ωi_m || ω_np_prime < self.ωi_0
            continue
        end

        A_pc = self.A_array[Ai_pc] 
        A_nc = self.A_array[Ai_nc] 

        δ_A = (δωi_c * (A_nc - A_pc) / (2*δωi_c + self.ωi_array[Ai_pc] - (-self.ωi_array[Ai_nc])))  #To conserve first moment

        A_pc_prime = A_pc + δ_A
        A_nc_prime = A_nc - δ_A

        A_pp = self.A_array[Ai_pp] 
        A_np = self.A_array[Ai_np] 

        δ_A = (δωi_p * (A_np - A_pp) / (2*δωi_p + self.ωi_array[Ai_pp] - (-self.ωi_array[Ai_np])))  #To conserve first moment

        A_pp_prime = A_pp + δ_A
        A_np_prime = A_np - δ_A

        if A_pc_prime < 0 || A_nc_prime < 0 || A_pp_prime < 0 || A_np_prime < 0
            continue
        end

        
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j]
                    + (A_pc_prime * self.Kp_array[j, ω_pc_prime + 1])
                    + (A_nc_prime * self.Kn_array[j, ω_nc_prime + 1])
                    - (A_pc * self.Kp_array[j, ω_pc + 1]) 
                    - (A_nc * self.Kn_array[j, ω_nc + 1])
                    + (A_pp_prime * self.Kp_array[j, ω_pp_prime + 1])
                    + (A_np_prime * self.Kn_array[j, ω_np_prime + 1])
                    - (A_pp * self.Kp_array[j, ω_pp + 1]) 
                    - (A_np * self.Kn_array[j, ω_np + 1]))
        end

        χ2_new = calc_χ2_new(self)
        
        P = exp((self.χ2 - χ2_new)/(2*θ))


        if rand() <= P
            self.ωi_array[Ai_pc] = ω_pc_prime
            self.ωi_array[Ai_nc] = ω_nc_prime
            self.A_array[Ai_pc] = A_pc_prime
            self.A_array[Ai_nc] = A_nc_prime

            self.ωi_array[Ai_pp] = ω_pp_prime
            self.ωi_array[Ai_np] = ω_np_prime
            self.A_array[Ai_pp] = A_pp_prime
            self.A_array[Ai_np] = A_np_prime

            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new

            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1.
        end
    end
    
    return accept_rate / (self.N_ω ÷ 2)
end



##########################################################################################
## Functions for performing above sampling routines and recording results  
##########################################################################################


# Function to run set of updates at a given θ
function run_updates(self::SAC, θ::Float64, transfer::Bool=true)

    sample_neg = true
    # Whether to sample pos. or pos. and neg. axes
    if self.symm == 1 || self.kernel_type == :bosonic
        pns = [1]
        sample_neg = false
        transfer = false
    else
        pns = [1, 2]
    end
   
    if transfer

        if self.indiv_update
            #a_r 6 shift weight between pos. and neg. continua
            ar = cont_Aω_transfer(self, θ)
            self.accept_rates[6] += ar

            #a_r 7 shift weight between pos. and neg. peaks
            ar = peak_Aω_transfer(self, θ)
            self.accept_rates[7] += ar
        end

        #a_r 1 combined update that shifts weight between pos. and neg. peaks and continua at same time
        ar = Aω_transfer(self, θ)
        self.accept_rates[11] += ar
       

    end


    for pn in pns
        ar = single_ω_move(self, θ, pn) #a_r 1 and 8 (pos and neg)
        self.accept_rates[1 + (pn-1)*7] += ar
        
        ar = double_ω_move(self, θ, pn) #a_r 2 and 9
        self.accept_rates[2 + (pn-1)*7] += ar
        
        ar = triple_ω_move(self, θ, pn) #a_r 3 and 10
        self.accept_rates[3 + (pn-1)*7] += ar

    end

    ar = single_ω_pmove(self, θ, 1) #a_r 4 move positions of single macro delta (+ axis)
    self.accept_rates[4] += ar
    
    if sample_neg
        ar = single_ω_pmove(self, θ, 2) #a_r 5 move positions of single macro delta (- axis)
        self.accept_rates[5] += ar
    end

    # self.ωi_array[self.peak_p] .= self.ωi_array[self.peak_n]
end

# Adjust update windows pased on acceptance rates 
# Windows are made larger/smaller to maintain ~50% acceptance rate for each update
function adjust_windows(self::SAC, steps::Int64, θ::Float64)
    
   
    for j=1:10
        
        self.accept_rates .= 0

        for i=1:(steps ÷ 10)

            self.Gbar = calc_Gbar(self)
            self.χ2 = calc_χ2(self)
            self.ωi_pp = peak_bound(self, 1)
            self.ωi_np = peak_bound(self, 2)
            run_updates(self, θ)
        end
        
        self.accept_rates /= (steps ÷ 10)

        for k=1:11
            if self.accept_rates[k] > 0.8
                self.update_windows[k] *= 2
            elseif self.accept_rates[k] < 0.2
                self.update_windows[k] /= 2
            elseif self.accept_rates[k] > 0.55
                self.update_windows[k] *= 1.2
            elseif self.accept_rates[k] < 0.45
                self.update_windows[k] /= 1.2
            end   
        end
        
    end


end

# Run sampling at fixed θ
# Can run as a single sweep in a bin for measuring fluctuations in edge location
function sample(self::SAC, steps::Int64, θ::Float64, bin::Bool = false)

    # if not a bin measurement, then initialize sampled spectra
    if bin == false 
        self.sampled_pspec .= 0
        self.sampled_cspec .= 0 
        transfer = true
    else
        transfer = false # Don't transfer weight between +/- if final sampling
    end

    self.sampled_χ2 = 0
    self.accept_rates .= 0
    self.edge_res[1, :] .= 0

    self.ωi_pp = peak_bound(self, 1)
    self.ωi_np = peak_bound(self, 2)
    
    for j=1:steps

        self.Gbar = calc_Gbar(self)
        self.χ2 = calc_χ2(self)

 
        run_updates(self, θ, transfer)
        
        for i in self.peak_p
            self.sampled_pspec[self.ωi_array[i]+1, 1] += self.A_array[i]
        end

        for i in self.peak_n
            self.sampled_pspec[self.ωi_array[i]+1, 2] += self.A_array[i]
        end

        for i in self.cont_p
            self.sampled_cspec[self.ωi_array[i]+1, 1] += self.A_array[i]
        end

        for i in self.cont_n
            self.sampled_cspec[self.ωi_array[i]+1, 2] += self.A_array[i]
        end

        self.ωi_pp = peak_bound(self, 1)
        self.ωi_np = peak_bound(self, 2)

        self.sampled_χ2 += self.χ2 # collected χ2 sample
        self.edge_res[1, 1] += self.ωi_pp # collect positive axis edge sample
        self.edge_res[1, 2] += self.ωi_np # collect positive axis edge sample
    end

    # if not a bin measurement, average accumulated spectra over number of sweeps
    if bin == false 
        self.sampled_pspec ./= steps
        self.sampled_cspec ./= steps
    end
    self.sampled_χ2 /= steps
    self.accept_rates ./= steps

    # For calculating error bar
    self.edge_res[1, :] ./= steps
    self.edge_res[2, :] .+= self.edge_res[1, :]
    self.edge_res[3, :] .+= self.edge_res[1, :] .^ 2

    if self.symm == 1 || self.kernel_type == :bosonic
        self.edge_res[:, 2] = self.edge_res[:, 1]
    end
end



function write_spec(self, n)

    spec_outfile = string(self.output_folder, "/sw", string(n, pad=3), ".csv")
    
    # Multiply back normalization factors
    f = 1. #additional normalzation factor for bosonic spectra
    for i=0:self.ωi_m 
        ω = self.δω * i
        if self.kernel_type == :bosonic
            f = 1 + exp(-self.β * ω)
        end
        self.sampled_pspec[i+1, :] *= (self.norm * pi)/f
        self.sampled_cspec[i+1, :] *= (self.norm * pi)/f
    end
    
    δω_conversion = Int64(round(self.δω_h/self.δω))
    N_h = Int64(self.ωi_m / δω_conversion)

    S_hist = 0.
    for pn in [1, 2]
        for i=0:(N_h-1)
            S_hist = sum(self.sampled_pspec[i*δω_conversion + 1:(i+1) * δω_conversion, pn])/self.δω_h
            self.sampled_pspec[i+1, pn] = S_hist

            S_hist = sum(self.sampled_cspec[i*δω_conversion + 1:(i+1) * δω_conversion, pn])/self.δω_h
            self.sampled_cspec[i+1, pn] = S_hist
        end
    end
    
    j0 = 1
    jf = N_h

    for i=N_h-1:-1:0
        jf=i
        if self.sampled_cspec[i+1, 1] > 1e-10
            break
        end
    end

    k0 = 1
    kf = N_h
    for i=N_h-1:-1:0
        kf=i
        if self.sampled_cspec[i+1, 2] > 1e-10
            break
        end
    end
    
     
    ω_hist = 0.
    S_p = 0.
    S_c = 0.
    f = open(spec_outfile, "w")
        println(f, "omega,S,S_p,S_c")

        # if bosonic spectrum, write S(-ω) = S(ω) exp(-βω)
        if self.kernel_type == :bosonic
            for i=jf:-1:j0
                ω_hist = self.δω_h * (i-1)
                S_p, S_c = self.sampled_pspec[i, 1], self.sampled_cspec[i, 1]
                S_hist = S_p + S_c

                println(f, -ω_hist, ",", S_hist*exp(-self.β*ω_hist), ",", S_p*exp(-self.β*ω_hist), ",", S_c*exp(-self.β*ω_hist))
            end
        # if symm, write S(-ω) = S(ω)
        elseif self.symm == 1
            for i=jf:-1:j0
                ω_hist = self.δω_h * (i-1)
                S_p, S_c = self.sampled_pspec[i, 1], self.sampled_cspec[i, 1]
                S_hist = S_p + S_c

                println(f, -ω_hist, ",", S_hist, ",", S_p, ",", S_c)
            end
        # othewise write sampled neg. spectrum 
        else
            for i=kf:-1:k0
                ω_hist = self.δω_h * (i-1)
                S_p, S_c = self.sampled_pspec[i, 2], self.sampled_cspec[i, 2]
                S_hist = S_p + S_c

                println(f, -ω_hist, ",", S_hist, ",", S_p, ",", S_c)
            end
        end

        # write sampled pos. spectrum 
        for i=j0:jf
            ω_hist = self.δω_h * (i-1)
            S_p, S_c = self.sampled_pspec[i, 1], self.sampled_cspec[i, 1]
            S_hist = S_p + S_c

            println(f, ω_hist, ",", S_hist, ",", S_p, ",", S_c)
        end

    close(f)

end

###########################################################################################
## Annealing Routines
###########################################################################################

# Higher temp intialization steps
function initial_sampling(self, θ_0)

    adjust_windows(self, self.anneal_steps, 10*θ_0)
    adjust_windows(self, self.anneal_steps, 5*θ_0)
    adjust_windows(self, self.anneal_steps, 2*θ_0)


end

# Run anneal, N_anneal temperature steps or until χ2_min converged 
# Can write spectrum each annealing step (wr = true) or only once after identifying optimal θ
function run_anneal(self, θ_0, wr=false)

    θ = θ_0

    for i=1:self.N_anneal
        
        # steps = floor(Int64, self.anneal_steps * sqrt(i)) # Ramp up number of steps per Θ
        steps = self.anneal_steps

        adjust_windows(self, steps, θ)
        sample(self, steps, θ)

        edge_p = self.ωi_pp * self.δω
        A0_p = sum(self.A_array[self.peak_p])
        Ac_p = sum(self.A_array[self.cont_p])
        if self.symm == 1 || self.kernel_type == :bosonic
            edge_n = -edge_p
            A0_n = A0_p
            Ac_n = Ac_p
        else
            edge_n = -self.ωi_np * self.δω
            A0_n = sum(self.A_array[self.peak_n])
            Ac_n = sum(self.A_array[self.cont_n])
        end
        
        
        open(self.anneal_file, "a") do f
            writedlm(f, [[i, map(x -> round(x, digits=8), 
                          [θ, self.χ2_min/self.N_τ, self.sampled_χ2/self.N_τ,
                           edge_p, edge_n, A0_p, A0_n, Ac_p, Ac_n])...]], ',')
        end

        # Write acceptance rates
            open(self.accept_rate_file, "a") do f
            output = cat([i], 
                         map(x -> round(x, digits=8), self.accept_rates[[1, 2, 3, 8, 9, 10]]),
                         map(x -> round(x, digits=8), self.update_windows[[1, 2, 8, 9]] .* self.δω),
                         map(x -> round(x, digits=8), self.accept_rates[[4, 5]]),
                         map(x -> round(x, digits=8), self.update_windows[[4, 5]] .* self.δω),
                         map(x -> round(x, digits=8), self.accept_rates[[6, 7, 11]]),
                         map(x -> round(x, digits=8), self.update_windows[[6, 7]].* self.δω),
                         [round(self.update_windows[11]* self.δω, digits=8)],
                         dims=1)
            writedlm(f, [output], ',')
        end
    

       
        if wr
            write_spec(self, i)
        end


        self.χ2_anneal[i] = self.sampled_χ2

        if (self.sampled_χ2 - self.χ2_min) < (self.tol * self.N_τ) 
            return
        else
            θ /= self.f_anneal
        end  

    end
end

# After identifying optimal θ from main anneal, run a final anneal from 10*θ_opt to θ_opt
# spectrum coresponding to a_criterion written to file
function final_anneal(self::SAC, θ_opt::Float64)

    write_log(self, "Beginning Fast Anneal to θ_opt.")
    for i=1:10
        θ = θ_opt * (11-i)
        steps = self.anneal_steps * i
        adjust_windows(self, steps, θ)
        sample(self, steps, θ)
    end
    write_log(self, "Fast Anneal to θ_opt Finished.")


    write_log(self, "Beginning Final Sampling at θ_opt.")
    # Accumulate ten bins of sampling sweeps for error bar
    N_bins = 10
    
    self.edge_res .= 0
    self.sampled_pspec .= 0
    self.sampled_cspec .= 0

    θ = θ_opt
    for n=1:N_bins
        sample(self, self.sample_steps, θ, true)
    end

    # average over all bins
    self.sampled_pspec ./= (N_bins * self.sample_steps)
    self.sampled_cspec ./= (N_bins * self.sample_steps)
    self.edge_res ./= N_bins
    self.edge_res[3, :] .= sqrt.(abs.(self.edge_res[3, :] .- self.edge_res[2, :] .^ 2)) ./ sqrt(N_bins - 1)
    self.edge_res .*= self.δω

    a = (self.sampled_χ2 - self.χ2_min) / sqrt(2 * self.χ2_min)
   
    A0_p = sum(self.A_array[self.peak_p])
    Ac_p = sum(self.A_array[self.cont_p])
    if self.symm == 1 || self.kernel_type == :bosonic
        A0_n = A0_p
        Ac_n = Ac_p
    else
        A0_n = sum(self.A_array[self.peak_n])
        Ac_n = sum(self.A_array[self.cont_n])
    end

    

    open(self.sample_file, "a") do f
        writedlm(f, [[0, 
                      round(a, digits=8), round(θ, digits=8), 
                      round(self.χ2_min/self.N_τ, digits=8), round(self.sampled_χ2/self.N_τ, digits=8),
                      round(self.edge_res[2, 1], digits=8), round(-self.edge_res[2, 2], digits=8), 
                      round(self.edge_res[3, 1], digits=8), round(self.edge_res[3, 2], digits=8), 
                      round(A0_p, digits=8), round(A0_n, digits=8),
                      round(Ac_p, digits=8), round(Ac_n, digits=8)]], ',')
    end

    open(self.accept_rate_file, "a") do f
        output = cat([0], 
                     map(x -> round(x, digits=8), self.accept_rates[[1, 2, 3, 8, 9, 10]]),
                     map(x -> round(x, digits=8), self.update_windows[[1, 2, 8, 9]] .* self.δω),
                     map(x -> round(x, digits=8), self.accept_rates[[4, 5]]),
                     map(x -> round(x, digits=8), self.update_windows[[4, 5]] .* self.δω),
                     map(x -> round(x, digits=8), self.accept_rates[[6, 7, 11]]),
                     map(x -> round(x, digits=8), self.update_windows[[6, 7]].* self.δω),
                     [round(self.update_windows[11], digits=8)],
                     dims=1)
        writedlm(f, [output], ',')
    end


    write_spec(self, 0)
    write_log(self, "Final Sampling at θ_opt Finished.")

end


function write_log(self, message, w = "a")
    f = open(self.output_folder * "/log.txt", w)
        println(f, string(now(), " - ", message))
    close(f)
end


# Reads all inputs from file "sac_peak.in," but also takes
# in A_0 and N_p to aid in scanning procedure
function run(A_0_in=false, N_p_in=false)

    in_file = readdlm("in_peak.in")

    N_ω, N_p = in_file[1, :]
    A_0, ω_m, δω, δω_h = in_file[2, :]
    θ_0, f_anneal, a_criterion = in_file[3, :]
    N_anneal, anneal_steps, sample_steps = in_file[4, :]
    G_file, output_folder = in_file[5, :]
    symm, kernel_type = in_file[6, 1], Symbol(in_file[6, 2])
    
   
    

    if A_0_in != false
        A_0 = A_0_in
    end
    if N_p_in != false
        N_p = N_p_in
    end

    if kernel_type == :bosonic
        symm = 0

    elseif symm == 1
        output_folder *= "_symm"
    end

    output_folder = output_folder * @sprintf("/Np_%02i/A0_%.3f", N_p, A_0)


    mkpath(output_folder)

    cp("in_peak.in", output_folder * "/in_peak.in", force = true)
    cp(G_file, output_folder * "/t.in", force = true)


    sac = SAC(N_ω=N_ω, N_p = N_p, ω_m=ω_m, δω=δω, δω_h=δω_h, A_0=A_0,
            anneal_steps=anneal_steps, sample_steps=sample_steps,
            θ_0=θ_0, N_anneal=N_anneal, f_anneal=f_anneal, a_criterion=a_criterion,
            G_file=G_file,
            anneal_file=string(output_folder, "/anneal.csv"),
            sample_file=string(output_folder, "/sample.csv"),
            accept_rate_file=string(output_folder, "/accept_rate.csv"), 
            output_folder=output_folder,
            kernel_type=kernel_type, symm=symm)


    # Initialization
    init_outfiles(sac)
    
    #STEP 1
    write_log(sac, "Beginning Initialization.", "w")
    
    read_G!(sac)
    initialize!(sac)

    sac.Gbar = calc_Gbar(sac)
    sac.Gbar_new = similar(sac.Gbar)
    sac.χ2 = calc_χ2(sac)
    sac.χ2_min = sac.χ2
    sac.update_windows .= sac.ω_window / 10
  
    θ = sac.θ_0



    
  
    write_log(sac, "Initialization Finished.")
    
    #STEP 2
   
  
    write_log(sac, "Beginning Initial Sampling.")
    initial_sampling(sac, θ)
    if sac.χ2_min > 1000 * sac.N_τ # Run initialization again if χ2 is too high first time 
        sac.indiv_update = true
        read_G!(sac)
        initialize!(sac)
        sac.Gbar = calc_Gbar(sac)
        sac.Gbar_new = similar(sac.Gbar)
        sac.χ2 = calc_χ2(sac)
        sac.χ2_min = sac.χ2
        sac.update_windows .= sac.ω_window / 10
        initial_sampling(sac, θ)
        # sac.indiv_update = false
    end
    write_log(sac, "Initial Sampling Finished.")

    #STEP 3
    write_log(sac, "Beginning Anneal.")
    run_anneal(sac, θ, false)
    write_log(sac, "Anneal Finished.")

    


    a_vals = (sac.χ2_anneal .- sac.χ2_min) ./ sqrt(2*sac.χ2_min)
    θ_vals = sac.θ_0 ./ (sac.f_anneal .^ collect(0:sac.N_anneal-1))

    θ_opt = θ_vals[argmin(abs.(a_vals .- sac.a_criterion))]

   
    f = open(string(output_folder, "/a_vals.csv"), "w")
        writedlm(f, [["a" "theta"]], ',')
        writedlm(f, [a_vals θ_vals], ',')
    close(f)
    

    sac.update_windows .= sac.ω_window / 10


    
    #STEP 4
    write_log(sac, "Beginning Final Anneal.")
    final_anneal(sac, θ_opt)
    write_log(sac, "Final Anneal Finished.")  
    
end


if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) > 0
        A_0, N_p = parse.(Float64, ARGS)
        run(A_0, N_p)
    else
        run()
    end
end


