# SAC code for FERMIONIC spectral function
# "Free/Unconstrained" Parameterization

using Dates, DelimitedFiles, LinearAlgebra, Printf

# Structure used to contain all the data for the routine
module SACs
Base.@kwdef mutable struct SAC
    par::Int64 # Delta function grid parameterization:
               # 1 = All deltas have equal amplitide d
               # 2 = Amplitudes of deltas are updated
               # 3 = Amplitudes of deltas are not all equal,
               #     but no ampltiude updates are carried out

    N_ω::Int64 # Number of delta functions
    ω_0::Float64 # Lower bound on deltas
    ω_m::Float64 # Upper bound on deltas
    ωi_0::Int64 = 0 # Index of lower bound on deltas
    ωi_m::Int64 = 0 # Index of upper bound on deltas
    δω::Float64 # Spacing of frequency grid for sampling 
    δω_h::Float64 # Spacing of frequency grid for output histogram of spectrum

    anneal_steps::Int64 # Number of sweeps preformed in each temperture step of the annealing processs
    sample_steps::Int64 # Number of sweeps preformed in the final sampling stage, usually larger than anneal_steps
    θ_0::Float64 # Initial sampling temperature
    N_anneal::Int64 # Number of temperature steps in anneal
    f_anneal::Float64 # Factor θ is divided by at each annealing steps: θ_{i+1} = θ_{i}/f_anneal
    f_final::Float64 # Factor θ is divided by at each annealing step in final anneal, before final sampling stage
    
    tol::Float64 = 1e-3 # Used to end annealing process early, is |χ^2_{i+1} - χ^2_{i}| < tol

    a1::Float64 # a value in optimal theta criterion, lower limit
    a2::Float64 # a value in optimal theta criterion, upper limit 

    G_file::String # File containing G(tau) data and covariance matrix

    β::Float64 = 0. # Inverse temperature
    N_τ::Int64 = 0 # Number of τ points
    norm::Float64 = 0. # Normalization of spectral function
    τ::Array{Float64, 1} = Array{Float64}(undef, 0) # Array of τ points
    G::Array{Float64, 1} = Array{Float64}(undef, 0) # Average of G(τ) bins
    σ_inv::Array{Float64, 1} = Array{Float64}(undef, 0) # Inverse of covariance matrix eigenvalues
    cov::Array{Float64, 2} = Array{Float64}(undef, 0, 0) # Covariance matrix
    ω_window::Float64 = 0. # Rough frequency scale deduced from G(τ)


    ωi_array::Array{Int64, 1} = Array{Int64}(undef, 0) # Delta function locations 
    A_array::Array{Float64, 1} = Array{Float64}(undef, 0) # Delta function amplitudes 
    K_array::Array{Float64, 2} = Array{Float64}(undef, 0, 0) # Kernel 

    Gbar::Array{Float64, 1} = Array{Float64}(undef, 0) #SAC G(τ) from ωi_array, A_array
    Gbar_new::Array{Float64, 1} = Array{Float64}(undef, 0) #SAC G(τ) after updatefrom ωi_array, A_array
    χ2::Float64 = 0. 
    χ2_new::Float64 = 0. # χ^2 after update
    χ2_min::Float64 = 0. # Minimum χ^2 during anneal
    
    χ2_anneal::Array{Float64, 1} = Array{Float64}(undef, 0) # χ^2 vs. Theta

    accept_rates::Array{Float64, 1} = zeros(6) # Acceptance rates for various updates
    update_windows::Array{Float64, 1} = zeros(4) # Frequency and ampltiude windows used for updates

    sampled_spec::Array{Float64, 1} = Array{Float64}(undef, 0) # Accumlated spectrum during sampling
    sampled_χ2::Float64 = 0. # Accumulated χ^2 during sampling

    anneal_file::String # Output file for annealing results
    accept_rate_file::String # Output file for acceptance rates
    sample_file::String # Output file for info on final sampling
    output_folder::String # Folder to write these files to


    kernel_type::Symbol # finiteT or zeroT or bosonic
        

    symm::Int64 # For fermionic spectral functions only!
                # Whether or not to impose particle-hole symettry for Hubbard like models
                # If set to true, A(-ω) = A(ω) and only the positive part of the spectrum is sampled
     
    
end
end
SAC = SACs.SAC




# Utility function used in updates
macro swap(x,y)
   quote
      local tmp = $(esc(x))
      $(esc(x)) = $(esc(y))
      $(esc(y)) = tmp
    end
end






# Read the input file "G_file"
# This file should be prepared by "make_tin.jl"
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
    
    ω_window = log(1/G[N_τ ÷ 2])/τ[N_τ ÷ 2]

    
    G_D = transpose(cov) * G # Transform into eigenbasis of covariance matrix
    
    self.τ = τ
    self.G = G_D
    self.σ_inv = σ_inv
    self.cov = cov
    self.ω_window = ω_window
    
end

#Initializing ωi_array, A_array, Kernel, etc.
function initialize!(self)

    if self.symm == 1
        self.ω_0 = 0
    end
    
    ωi_0 = Int64(floor(self.ω_0/self.δω)) # Converting ω_0 to integer
    ωi_m = Int64(ceil(self.ω_m/self.δω))
    
    ω_window = self.ω_window/self.δω
    
    ωi_array = zeros(Int64, self.N_ω)
    A_array = zeros(Float64, self.N_ω)

    if ωi_0 < 0
        N_ωi = Int64(ωi_m-ωi_0+1) # Slots for negative and positive frquencies
        ω_start = ωi_0
        ∆ = ((ωi_m - (ωi_0))) ÷ (self.N_ω) # Spacing for initial distribution of deltas on frequency grid
    else
        N_ωi = Int64(ωi_m+1) 
        ω_start = 0
        ∆ = (ωi_m - (ωi_0)) ÷ (self.N_ω)
    end

    self.sampled_spec = zeros(Float64, N_ωi)

    self.χ2_anneal = zeros(self.N_anneal)
    
    

    for i=1:self.N_ω
   
        ωi_array[i] = ω_start + (∆*i)
        
        if self.par == 1 || self.par == 2 # If par = 1,2 then all deltas are intialized with same ampltiude
            A_array[i] = 1
        elseif self.par == 3 # If par = 3 then deltas are intialized with increasing ampltiude
             A_array[i] = i
        end
    end


    if self.symm == 1
        A_array ./= (2*sum(A_array)) # If symm, then only positive side of frequency axis is sampled,
                                     # and is copied to negative frequency axis.
                                     # Need to divide amplitude by two to maintain normaliztion
                                     # Σ A_i = 1
    else
        A_array ./= (sum(A_array)) # Normalization Σ A_i = 1
    end

    
    #Initializing Fermionic Kernel

    if self.kernel_type == :finiteT
        K_function = finiteT_K
    elseif self.kernel_type == :zeroT
        K_function = zeroT_K
    elseif self.kernel_type == :bosonic
        K_function = bosonic_K
    else
        throw("Invalid Kernel type.")
    end

    
    K_array = zeros(Float64, (self.N_τ, N_ωi))
    K_vals = zeros(Float64, self.N_τ)
    for i=ωi_0:ωi_m
        ω = ((i + 0.5) * self.δω)
        K_vals .= K_function.(ω, self.τ, self.β)

        #If symm, then kernel is adjusted to account for S(-ω) = S(ω) 
        if self.symm == 1  
            K_vals .+=  K_function.(-ω, self.τ, self.β)
        end
      
       

        K_array[:, i-ωi_0+1] .= K_vals
        K_array[:, i-ωi_0+1] .= transpose(self.cov) * K_array[:, i-ωi_0+1] # Transform into eigenbasis of covariance matrix
    end
    
    self.ωi_array, self.A_array, self.K_array = ωi_array, A_array, K_array
    self.ω_window, self.ωi_0, self.ωi_m = ω_window, ωi_0, ωi_m
    
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


# G_bar = G_SAC(τ) = Σ K(ω_i, τ) A(ω_i)
function calc_Gbar(self::SAC)
    Gbar = self.K_array[:, self.ωi_array .+ 1 .- self.ωi_0] * self.A_array
    return Gbar 
end

# χ^2 = Σ (G_bar - G) / σ^2
function calc_χ2(self::SAC)
    χ2 = 0
    @inbounds @simd for i=1:self.N_τ
        χ2 += ((self.Gbar[i] - self.G[i]) * self.σ_inv[i]) ^ 2
    end
    return χ2
end

function calc_χ2_new(self::SAC)
    χ2 = 0
    @inbounds @simd for i=1:self.N_τ
        χ2 += ((self.Gbar_new[i] - self.G[i]) * self.σ_inv[i]) ^ 2
    end
    return χ2
end

##########################################################################################
##########################################################################################



# Sampling Frequencies
##########################################################################################

# Adjust position of a single delta function
function single_ω_move(self::SAC, θ::Float64)
    
    accept_rate = 0
    
    update_window = self.update_windows[1] # Magnitude of change in freuqnecy


    ω_1 = -1

    for i=1:self.N_ω
        ωi = rand(1:self.N_ω)
        δωi = 1 + Int64(floor(rand()*update_window))
        
        if rand() < 0.5
            ω_1 = self.ωi_array[ωi] + δωi
        else
            ω_1 = self.ωi_array[ωi] - δωi
        end
        
        if ω_1 < self.ωi_0 
            continue
        elseif ω_1 > self.ωi_m
            continue
        end
        
        A = self.A_array[ωi]
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
              (A * (self.K_array[j, ω_1 - self.ωi_0 + 1] - self.K_array[j, self.ωi_array[ωi] - self.ωi_0 + 1])))
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
    
    return accept_rate/self.N_ω # Return accept of full sweep
end


# Adjust position of a pair of delta functions
function double_ω_move(self::SAC, θ::Float64)

    accept_rate = 0
    update_window = self.update_windows[2] # Magnitude of change in pair of freuqnecies

    ω_1 = -1
    ω_2 = -1

    for i=1:(self.N_ω ÷ 2)
        ωi_1 = rand(1:self.N_ω)
        ωi_2 = ωi_1
        while ωi_2 == ωi_1
            ωi_2 =  rand(1:self.N_ω)
        end
        
        δωi = 1 + Int64(floor(rand()*update_window))
        
        if rand() < 0.5
            ω_1 = self.ωi_array[ωi_1] + δωi
            ω_2 = self.ωi_array[ωi_2] - δωi
        else
            ω_1 = self.ωi_array[ωi_1] - δωi
            ω_2 = self.ωi_array[ωi_2] + δωi
        end
        
        if ω_1 < self.ωi_0
            continue
        elseif ω_2 < self.ωi_0
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
                A1 * (self.K_array[j, ω_1 - self.ωi_0 + 1] - self.K_array[j, self.ωi_array[ωi_1] - self.ωi_0  + 1]) +
                A2 * (self.K_array[j, ω_2 - self.ωi_0 + 1] - self.K_array[j, self.ωi_array[ωi_2] - self.ωi_0  + 1]))
        end

        
        χ2_new  = calc_χ2_new(self)
        #println(χ2_new, " ", χ2)
        
        @fastmath P = exp((self.χ2 - χ2_new)/(2*θ))
        #println(P)
        
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
    return accept_rate/(self.N_ω ÷ 2) # Return accept of full sweep
    
    
    
end


# Adjust position of a trio of delta functions
function triple_ω_move(self::SAC, θ::Float64)
    #println("triple")
    
    accept_rate = 0

    ω_1 = -1 
    ω_2 = -1
    ω_3 = -1

    
    for i=1:(self.N_ω ÷ 3)
        ωi_1 = rand(1:self.N_ω)
        ωi_2 = ωi_1
        while ωi_2 == ωi_1
            ωi_2 =  rand(1:self.N_ω)
        end
        ωi_3 = ωi_2
        while ωi_3 == ωi_1 || ωi_3 == ωi_2
            ωi_3 =  rand(1:self.N_ω)
        end
        
        # No update window used, ∆ω set to conserve first moment

        δωi = (self.ωi_array[ωi_2] + self.ωi_array[ωi_3] - (2 * self.ωi_array[ωi_1])) ÷ 3
        
        ω_1 = self.ωi_array[ωi_1] + 2*δωi
        ω_2 = self.ωi_array[ωi_2] - δωi
        ω_3 = self.ωi_array[ωi_3] - δωi
        
        if ω_1 < self.ωi_0
            continue
        elseif ω_2 < self.ωi_0
            continue
        elseif ω_3 < self.ωi_0
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
                A1 * (self.K_array[j, ω_1 - self.ωi_0 + 1] - self.K_array[j, self.ωi_array[ωi_1] - self.ωi_0 + 1]) +
                A2 * (self.K_array[j, ω_2 - self.ωi_0 + 1] - self.K_array[j, self.ωi_array[ωi_2] - self.ωi_0 + 1]) +
                A3 * (self.K_array[j, ω_3 - self.ωi_0 + 1] - self.K_array[j, self.ωi_array[ωi_3] - self.ωi_0 + 1]))
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
    return accept_rate/(self.N_ω ÷ 3) # Return accept of full sweep
  
end



# Sampling Amplitudes
##########################################################################################

# Adjust position and frequency of pair of delta functions, to conserve first moment
function A_ω_move(self::SAC, θ::Float64)

    accept_rate = 0
   
    update_window = self.update_windows[4]

    ω_1 = -1 
    ω_2 = -1
    
    for i=1:(self.N_ω ÷ 2)
        ωi_1 = rand(1:self.N_ω)
        ωi_2 = ωi_1
        while ωi_2 == ωi_1
            ωi_2 =  rand(1:self.N_ω)
        end

        δωi = 1 + Int64(floor(rand()*update_window))

        if rand() < 0.5
            ω_1 = self.ωi_array[ωi_1] + δωi
            ω_2 = self.ωi_array[ωi_2] - δωi
        else
            ω_1 = self.ωi_array[ωi_1] - δωi
            ω_2 = self.ωi_array[ωi_2] + δωi
        end

        if ω_1 < self.ωi_0
            continue
        elseif ω_2 < self.ωi_0
            continue
        elseif ω_1 > self.ωi_m
            continue
        elseif ω_2 > self.ωi_m
            continue
        end

        A1 = self.A_array[ωi_1] 
        A2 = self.A_array[ωi_2] 

        δ_A = (δωi * (A1 - A2) / (2*δωi + self.ωi_array[ωi_1] - self.ωi_array[ωi_2]))  #To conserve first moment

    
        A1_prime = A1 + δ_A
        A2_prime = A2 - δ_A

        if A1_prime < 0
            continue
        elseif A2_prime < 0
            continue
        end

        
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
                    (A1_prime * self.K_array[j, ω_1 - self.ωi_0 + 1]) + 
                    (A2_prime * self.K_array[j, ω_2 - self.ωi_0 + 1]) -
                    (A1 * self.K_array[j, self.ωi_array[ωi_1] - self.ωi_0 + 1]) - 
                    (A2 * self.K_array[j, self.ωi_array[ωi_2] - self.ωi_0 + 1]))
        end

        χ2_new = calc_χ2_new(self)
        
        P = exp((self.χ2 - χ2_new)/(2*θ))


        if rand() <= P
            self.ωi_array[ωi_1] = ω_1
            self.ωi_array[ωi_2] = ω_2
            self.A_array[ωi_1] = A1_prime
            self.A_array[ωi_2] = A2_prime

            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new

            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 2
        end
    end
    if self.symm == 1
        self.A_array ./= 2*sum(self.A_array)
    else
        self.A_array ./= sum(self.A_array)
    end

    return accept_rate/(self.N_ω ÷ 2)
end


# Swap amplitudes of pair of delta functions
function A_swap(self::SAC, θ::Float64)

    accept_rate = 0
    Ai_1 = -1
    Ai_2 = -1

    

    for i=1:self.N_ω
        Ai_1 = rand(1:self.N_ω)
        Ai_2 = Ai_1
        while self.ωi_array[Ai_1] == self.ωi_array[Ai_2]
            Ai_2 =  rand(1:self.N_ω)
        end

        δA_1 = self.A_array[Ai_2] - self.A_array[Ai_1]
        δA_2 = self.A_array[Ai_1] - self.A_array[Ai_2]
        
       @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
                (δA_1 * self.K_array[j, self.ωi_array[Ai_1] - self.ωi_0 + 1]) + 
                (δA_2 * self.K_array[j, self.ωi_array[Ai_2] - self.ωi_0 + 1]))
        end
        
        χ2_new = calc_χ2_new(self)
        
        P = exp((self.χ2 - χ2_new)/(2*θ))


        if rand() <= P
            self.A_array[Ai_1] += δA_1
            self.A_array[Ai_2] += δA_2

            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new

            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1
        end
    end
    return accept_rate/(self.N_ω)
end

# Update ampltiude of "adjacent" delta functions
# Currently not implemented
function neighbor_A_move(self::SAC, θ::Float64)
    
    accept_rate = 0
    Ai_1 = -1
    Ai_2 = -1
    
    for i=1:self.N_ω-1
        Ai_1 = rand(1:self.N_ω-1)
        Ai_2 = Ai_1 + 1
        
        m0 = self.A_array[Ai_1] + self.A_array[Ai_2]
        
        r = rand()

        δA_1 = (r * m0) - self.A_array[Ai_1]
        δA_2 = ((1 - r) * m0) - self.A_array[Ai_2]


        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
                (δA_1 * self.K_array[j, self.ωi_array[Ai_1] - self.ωi_0 + 1]) + 
                (δA_2 * self.K_array[j, self.ωi_array[Ai_2] - self.ωi_0 + 1]))
        end
        
        χ2_new = calc_χ2_new(self)
        
        P = exp((self.χ2 - χ2_new)/(2*θ))


        if rand() <= P
            self.A_array[Ai_1] += δA_1
            self.A_array[Ai_2] += δA_2

            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new

            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1
        end
    end
    
    if self.symm == 1
        self.A_array ./= 2*sum(self.A_array)
    else
        self.A_array ./= sum(self.A_array)
    end

    return accept_rate/(self.N_ω-1)
end

# Update ampltiudes of any pair of delta functions
function double_A_move(self::SAC, θ::Float64)

    accept_rate = 0
    Ai_1 = -1
    Ai_2 = -1
    
    
    for i=1:self.N_ω
        Ai_1 = rand(1:self.N_ω)
        Ai_2 = Ai_1
        while self.ωi_array[Ai_1] == self.ωi_array[Ai_2]
            Ai_2 = rand(1:self.N_ω)
        end
        
        m0 = self.A_array[Ai_1] + self.A_array[Ai_2]
        
        r = rand()

        δA_1 = (r * m0) - self.A_array[Ai_1]
        δA_2 = ((1 - r) * m0) - self.A_array[Ai_2]
        
        
        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
                (δA_1 * self.K_array[j, self.ωi_array[Ai_1] - self.ωi_0 + 1]) + 
                (δA_2 * self.K_array[j, self.ωi_array[Ai_2] - self.ωi_0 + 1]))
        end

        
        χ2_new = calc_χ2_new(self)
        
        P = exp((self.χ2 - χ2_new)/(2*θ))

        if rand() <= P
            self.A_array[Ai_1] += δA_1
            self.A_array[Ai_2] += δA_2

            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new

            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1
        end
    end
    
    if self.symm == 1
        self.A_array ./= 2*sum(self.A_array)
    else
        self.A_array ./= sum(self.A_array)
    end

    return accept_rate/(self.N_ω)
end

# Update ampltiudes of any trio of delta functions
# There have been some issues with this update, currently not implemented
function triple_A_move(self::SAC, θ::Float64)

    accept_rate = 0
    Ai_1 = -1
    Ai_2 = -1
    Ai_3 = -1
    

    for i=1:(self.N_ω ÷ 2)
        Ai_1 = rand(1:self.N_ω)
        Ai_2 = Ai_1
        while self.ωi_array[Ai_1] == self.ωi_array[Ai_2]
            Ai_2 = rand(1:self.N_ω)
        end
        
        if self.ωi_array[Ai_2] < self.ωi_array[Ai_1]
            @swap(Ai_1, Ai_2)
        end
        
        Ai_3 = Ai_2
        while self.ωi_array[Ai_3] == self.ωi_array[Ai_1] || self.ωi_array[Ai_3] == self.ωi_array[Ai_2] 
            Ai_3 = rand(1:self.N_ω)
        end
        
        if self.ωi_array[Ai_3] < self.ωi_array[Ai_1]
            @swap(Ai_3, Ai_2)
            @swap(Ai_2, Ai_1)
        elseif self.ωi_array[Ai_3] < self.ωi_array[Ai_2]
            @swap(Ai_3, Ai_2)
        end
        
        ω_1 = self.ωi_array[Ai_1] + 0.5 - self.ωi_0 
        ω_2 = self.ωi_array[Ai_2] + 0.5 - self.ωi_0 
        ω_3 = self.ωi_array[Ai_3] + 0.5 - self.ωi_0 

        A1 = self.A_array[Ai_1]
        A2 = self.A_array[Ai_2]
        A3 = self.A_array[Ai_3]
        
        m0 = A1 + A2 + A3
        m1 = ((A1 * ω_1) + (A2 * ω_2) + (A3 * ω_3))
        
        a1 = max(0, ((m0 * ω_2) - m1)/(ω_2 - ω_1))
        a2 = min(m0, ((m0 * ω_3) - m1)/(ω_3 - ω_1))
        
        r = rand()

        A1_prime = a1 + (r * (a2 - a1))
        A2_prime = ((m0 * ω_3) - m1 - (A1_prime * (ω_3 - ω_1))) / (ω_3 - ω_2)
        A3_prime = m0 - A1_prime - A2_prime
        

        @inbounds @simd for j=1:self.N_τ
            @views self.Gbar_new[j] = (self.Gbar[j] +
                ((A1_prime - A1)  * self.K_array[j, self.ωi_array[Ai_1] - self.ωi_0 + 1]) + 
                ((A2_prime - A2)  * self.K_array[j, self.ωi_array[Ai_2] - self.ωi_0 + 1]) +
                ((A3_prime - A3)  * self.K_array[j, self.ωi_array[Ai_3] - self.ωi_0 + 1]))
        end
        
        χ2_new = calc_χ2_new(self)
        
        P = exp((self.χ2 - χ2_new)/(2*θ))

        if rand() <= P
            self.A_array[Ai_1] += A1_prime
            self.A_array[Ai_2] += A2_prime
            self.A_array[Ai_3] += A3_prime

            self.Gbar .= self.Gbar_new
            self.χ2 = χ2_new

            if self.χ2 < self.χ2_min
                self.χ2_min = χ2_new
            end
            accept_rate += 1
        end
    end
    
    if self.symm == 1
        self.A_array ./= 2*sum(self.A_array)
    else
        self.A_array ./= sum(self.A_array)
    end
    
    return accept_rate/(self.N_ω÷ 2)
end




##########################################################################################
###########################################################################################

function run_updates(self::SAC , θ::Float64)
    
    # If par = 1,3 only preform frequency updates
    if self.par == 1 || self.par == 3
        ar = single_ω_move(self, θ) #a_r 1
        self.accept_rates[1] += ar
        
        ar = double_ω_move(self, θ) #a_r 2
        self.accept_rates[2] += ar
        
        ar = triple_ω_move(self, θ) #a_r 3
        self.accept_rates[3] += ar
    
    # If par = 2 preform both frequency and amplitude updates  
    elseif self.par == 2
        ar = single_ω_move(self, θ) #a_r 1
        self.accept_rates[1] += ar
        
        ar = double_ω_move(self, θ) #a_r 2
        self.accept_rates[2] += ar
        
        ar = triple_ω_move(self, θ) #a_r 3
        self.accept_rates[3] += ar



        ar = A_ω_move(self, θ) #a_r 4
        self.accept_rates[4] += ar

        ar = double_A_move(self, θ) #a_r 5
        self.accept_rates[5] += ar

        #ar = triple_A_move(self, θ) #a_r 6
        #self.accept_rates[6] += ar
    end

end


# Adjust update windows pased on acceptance rates 
# Windows are made larger/smaller to maintain ~50% acceptance rate for each update
function adjust_windows(self::SAC, steps::Int64, θ::Float64)
    
    # Preform some 10 sets of sweeps sweeps to accumulate acceptance rates
    for j=1:10
        
        self.accept_rates .= 0

        for i = 1:(steps ÷ 10)


            Gbar = calc_Gbar(self)
            χ2 = calc_χ2(self)

            run_updates(self, θ)
        end
        

        self.accept_rates ./= (steps ÷ 10)

        for k=1:4
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
function sample(self::SAC, steps::Int64, θ::Float64)
    

    self.sampled_spec .= 0 
    self.sampled_χ2 = 0
    self.accept_rates .= 0


    for j=1:steps

        Gbar = calc_Gbar(self)
        χ2 = calc_χ2(self)
        
        run_updates(self, θ)
        
        self.sampled_spec[self.ωi_array .+ 1 .- self.ωi_0] .+= self.A_array # Accumulate using updated ωi_array
        self.sampled_χ2 += self.χ2 # Collected χ2 sample
    end

    
    # Average over total number of update sweeps
    self.sampled_spec ./= steps
    self.sampled_χ2 /= steps
    self.accept_rates ./= steps
end



function init_outfiles(self)

    #All outputted as CSV's
   
    f = open(self.anneal_file, "w")
    writedlm(f, [["i", "theta", "chi2_min", "chi2_avg"]], ',')
    close(f)

    # Accept rates are: single ω, double ω, triple ω, A+ω, double A, triple A
    # omega ranges are for single ω, double ω, and A+ω updates
    f = open(self.accept_rate_file, "w")
    writedlm(f, [["i", "ar_1", "ar_2", "ar_3", "ar_4", "ar_5", "ar_6", "omega_range_1", "omega_range_2", "omega_range_4" ]], ',')
    close(f)


    f = open(self.sample_file, "w")
    writedlm(f, [["i", "a", "theta", "chi2_min", "chi2_avg"]], ',')
    close(f)


    f = open(string(self.output_folder, "/log.txt"), "w")
    close(f)


end


# Write rebinned histogram of accumulated spectrum from sampling
function write_spec(self, n)
    
    
    spec_outfile = string(self.output_folder, "/sw", string(n, pad=3), ".csv")
    
    δω_conversion = Int64(round(self.δω_h/self.δω))
    N_h = Int64((self.ωi_m - self.ωi_0) / δω_conversion)
    sampled_spec_hist = zeros(Float64, N_h)
    
    # Multiply back normalization factors
    f = 1. #additional normalzation factor for bosonic spectra
    for i=self.ωi_0:self.ωi_m
        ω = self.δω * i
        if self.kernel_type == :bosonic
            f = 1 + exp(-self.β * ω)
        end
        self.sampled_spec[i-self.ωi_0 + 1] *= (self.norm * pi)/f
    end

    # Rebin accumlated histogram 
    S_hist = 0.
    for i=1:N_h
        S_hist = sum(self.sampled_spec[1+(i-1)*δω_conversion:1+(i)*δω_conversion])/self.δω_h
        self.sampled_spec[i] = S_hist
    end

    j0 = 1
    jf = N_h

    for i=N_h-1:-1:0
        jf=i
        if self.sampled_spec[i] > 1e-10
            break
        end
    end


    # Write to file
    
    ω_hist = 0.
    f = open(spec_outfile, "w")
        println(f, "omega,S")

        # if bosonic spectrum, also write S(-ω) = S(ω) exp(-βω)
        if self.kernel_type == :bosonic
            for i=jf:-1:j0
                ω_hist = (self.δω_h * (i-1)) + (self.ωi_0*self.δω)
                println(f, -ω_hist, ",", self.sampled_spec[i]*exp(-self.β*ω_hist))
            end
        end

        for i=j0:jf

            ω_hist = (self.δω_h * (i-1)) + (self.ωi_0*self.δω)
            println(f, ω_hist, ",", self.sampled_spec[i])
        end
    close(f)

   
  
    
end


##########################################################################################
###########################################################################################



## a) High initial temperature, fast anneal to equillibrate, from θ=100*θ_0 to θ_0 in 10 steps
##    using 1/5 annealing steps per temp (currently implemented)
## b) Main anneal from θ_0 to θ_0 / (f_anneal)^N_anneal, storing mean χ2 and tracking χ^2 mean,
##    using constant number of steps θ (can also ramp up steps)
## c) Final anneal from 10*θ(a=a1) to θ(a=a1) using linear ramp of steps per θ,
##    after which sampled at fixed temperature to accumulate spectrum


# Fast anneal from very high temp to equillibrate delta functions after initialization 
function fast_anneal(self)

    for i=1:10
        θ = self.θ_0 * (11-i)^2
        
        adjust_windows(self, self.anneal_steps ÷ 2, θ)
        
        sample(self, self.anneal_steps ÷ 2, θ)

        # χ2 is normalized by number of τ points in output files
        open(self.anneal_file, "a") do f
           writedlm(f, [[i, θ, self.χ2_min/self.N_τ, self.sampled_χ2/self.N_τ]], ',')
        end
        
        open(self.accept_rate_file, "a") do f
            output = cat([i], self.accept_rates,self.update_windows[[1, 2, 4]].* self.δω, dims=1)
            writedlm(f, [output], ',')
        end
         
    end

end

# Run anneal, N_anneal temperature steps or until χ2_min converged 
function main_anneal(self)

    θ = self.θ_0

    for i=1:self.N_anneal
        #steps = Int64(floor(self.anneal_steps * sqrt(i))) # Ramp up number of steps per Θ
        steps = self.anneal_steps
        
        adjust_windows(self, steps, θ)
        
        sample(self, steps, θ)
        

        open(self.anneal_file, "a") do f
           writedlm(f, [[i, θ, self.χ2_min/self.N_τ, self.sampled_χ2/self.N_τ]], ',')
        end
        
        open(self.accept_rate_file, "a") do f
            output = cat([i], self.accept_rates,self.update_windows[[1, 2, 4]].* self.δω, dims=1)
            writedlm(f, [output], ',')
        end
        
        self.χ2_anneal[i] = self.sampled_χ2
        
        # Uncomment this line if you want to output the spectrum at each temperature in anneal
        # instead of just during final sampling step 
        #write_spec(self, i)
        
        # Exit main anneal is sufficiently close to χ2_min
        if (self.sampled_χ2 - self.χ2_min) < (self.tol * self.N_τ) 
            return
        else
            θ /= self.f_anneal
        end  

    end
end



# After identifying optimal θ from main anneal, run a final anneal from 10*θ_opt to θ_opt
# if a1 == a2, onlly one spectrum is output, otherwise spectra at a series of decreasing 
# a values are output

function final_anneal(self, θ_opt)

    for i=1:10
        θ = θ_opt * (11-i)

        steps = self.anneal_steps * i
        
        adjust_windows(self, steps, θ)

        sample(self, steps, θ)
         
    end

    if self.a1 == self.a2 # If a1 = a2 just run sampling once at θ_opt 
        N_final_anneal = 1
    else
        N_final_anneal = 20 # Max number of annealing steps between a1 and a2
                            # will lower θ from θ_opt until a = a1
    end
    
    θ = θ_opt
    for n=1:N_final_anneal

        sample(self, self.sample_steps, θ)
        
        # Calculate a value after sampling
        a = (self.sampled_χ2 - self.χ2_min) / sqrt(2 * self.χ2_min)

        
        open(self.sample_file, "a") do f
            writedlm(f, [[n-1, a, θ, self.χ2_min/self.N_τ, self.sampled_χ2/self.N_τ]], ',')
        end

        write_spec(self, n-1) # Output spectra with naming concention sw00n.csv
                              # where n increases from 0 (correpsonding to a2) until a1 is reached

        if a < self.a1
            return
        end

        θ /= self.f_final
    end

end




function write_log(self, message, w_a = "a")
    f = open(self.output_folder * "/log.txt", w_a)
        println(f, string(now(), " - ", message))
    close(f)
end





function run()

    # Input parameters for program
    in_file = readdlm("in_free.in")

    par = in_file[1,1]
    N_ω, ω_0, ω_m, δω, δω_h = in_file[2, :]
    θ_0, f_anneal, f_final, a1, a2 = in_file[3, :]
    N_anneal, anneal_steps, sample_steps = in_file[4, :]
    G_file, output_folder = in_file[5, :]
    symm, kernel_type = in_file[6, 1], Symbol(in_file[6, 2])
    
    if kernel_type == :bosonic
        ω_0 = 0.
        symm = 0
    elseif symm == 1
        output_folder *= "_symm"
        ω_0 = 0.
    end

    # Make output folder, in case path doesn't already exist (won't error if it does though)
    mkpath(output_folder)
    # Copy input parameters and G_file into output folder for reference
    cp("in_free.in", output_folder * "/in_free.in", force = true)
    cp(G_file, output_folder * "/t.in", force = true)



   
    # Construct struct using all provided info
    sac = SAC(par=par, N_ω=N_ω, ω_0=ω_0, ω_m =ω_m,δω=δω, δω_h=δω_h,
            anneal_steps=anneal_steps, sample_steps=sample_steps, θ_0=θ_0,
            N_anneal=N_anneal, f_anneal=f_anneal, f_final=f_final, a1=a1, a2=a2,
            G_file=G_file,
            anneal_file= output_folder * "/anneal.csv",
            accept_rate_file= output_folder * "/accept_rate.csv", 
            sample_file= output_folder *  "/sample.csv",
            output_folder=output_folder,
            symm=symm, kernel_type=kernel_type)


    ##########################################

    #STEP 1: Initialization

    write_log(sac, "Beginning Initialization.", "w")

    init_outfiles(sac) # Intialize output files
   
    read_G!(sac) # Read G(τ) and covariance matrix
    initialize!(sac) # Initialize spectrum

    sac.Gbar = calc_Gbar(sac)
    sac.Gbar_new = similar(sac.Gbar)
    sac.χ2 = calc_χ2(sac)
    sac.χ2_min = sac.χ2

    sac.update_windows .= sac.ω_window # Set default update window size
    
    write_log(sac, "Initialization Finished.")

    ##########################################

    #STEP 2: Equillibration

    write_log(sac, "Beginning Equillibration.")

    fast_anneal(sac)
    #equillibration(sac) # Not currently implemented

    write_log(sac, "Equillibration Finished.")

    ##########################################

    #STEP 3: Main anneal

    write_log(sac, "Beginning Main Anneal.")

    main_anneal(sac)

    write_log(sac, "Main Anneal Finished.")

    ##########################################
    
    #STEP 4: Final anneal

    write_log(sac, "Beginning Final Anneal.")

    a_vals = (sac.χ2_anneal .- sac.χ2_min) ./ sqrt(2*sac.χ2_min) # Calculate a values vs. θ to identify optimal θ
    θ_vals = sac.θ_0 ./ (sac.f_anneal .^ collect(0:sac.N_anneal-1))

    θ_opt = θ_vals[argmin(abs.(a_vals .- sac.a2))] # θ_opt is one where a is closest to a2

    sac.update_windows .= sac.ω_window # Reset update windows for final anneal
    

    final_anneal(sac, θ_opt)

    write_log(sac, "Final Anneal Finished.")

end


if abspath(PROGRAM_FILE) == @__FILE__
    run()
end


