using DelimitedFiles, LinearAlgebra, Statistics, Printf, Dates

module SACs
Base.@kwdef mutable struct SAC

   
    N_ω::Int64 # total number of δ's
    N_e::Int64 # number of δ's in edge
    A_c::Float64 = 0. # weight fraction in continuum
    A_r::Float64 # weight fraction in rightward decaying edge (1-A_r in leftward)
    d0::Float64 = 0. # min frequency gap at edge (ω_1 - ω_0) for quenching sharp edge
    ω_0::Float64 # min frequency
    ω_floor::Array{Float64, 1} = [0., 0.] # lowest frequency in sampling, used to fix edge location
    ω_m::Float64 # max frequency
    ωi_0::Int64 = 0 # min frequency index
    ωi_m::Int64 = 0 # max frequency index
    N_ωi::Int64 = 0 # total number of frequency indices = ωi_m - ωi_0
    ω_b0::Array{Float64, 1} = zeros(2) # frequency bounds of edges
    ω_bf::Array{Float64, 1} = zeros(2) #frequency bounds of tails
    δω::Float64 # frequency spacing in kernel grid and accumulated histograms
    δω_h::Float64 # frequency spacing in written histogram

    c::Array{Float64, 1} # exponent in power-law adjustment of weights, (1-c)/2 = p for edge diverging as (ω - ω0)^(-p)
    n0::Array{Float64, 1} = Array{Float64}(undef, 1) # border in amplitude form
    ε_0::Array{Float64, 1} = Array{Float64}(undef, 1) # rounding factor in amplitude form
    ∆n0::Array{Float64, 1} = [2., 2.] # step in changing n0
    ∆ε_0::Array{Float64, 1} = [0.1, 0.1] # step in changing ε_0

    anneal_steps::Int64 # base number of update sweeps per annealing step
    sample_steps::Int64 # number of updates sweeps in the final sampling stage
    bins::Int64 # number of bins used to determine std's of χ2, edge location, etc. in annealing
    θ_0::Float64 # initial ssampling temp
    N_anneal::Int64 # number of annealing steps 
    f_anneal::Float64 # theta ramp factor: θ → θ * f_ammeal 
    tol::Float64 = 1.0e-3 # stop anneal when  ∆χ2 < tol

    a_criterion::Float64 # a value to output spectra for in sampling stage
    
    G_file::String # name of input G(τ) data file

    # The following are read in from G_file
    β::Float64 = 0. # inverse temp of G(τ) data
    N_τ::Int64 = 0 # number of τ points
    norm::Float64 = 0. # G(0)
    τ::Array{Float64, 1} = Array{Float64}(undef, 1) # τ grid
    G::Array{Float64, 1} = Array{Float64}(undef, 1) # G(τ)
    σ_inv::Array{Float64, 1} = Array{Float64}(undef, 1) # covariance matric eigenvalues
    cov::Array{Float64, 2} = Array{Float64}(undef, 1, 1) # covariance matric eigenvectors
    edge_guess::Float64 = 0. # initial guess for edge location based on G(τ_max)

    ω_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1) # sampled frequencies
    d_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1) # distances between frequencies
    A_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1) # sampled amplitudes
    

    K_array::Array{Float64, 3} = Array{Float64}(undef, 1, 1, 1) # kernel
    dK_array::Array{Float64, 3} = Array{Float64}(undef, 1, 1, 1) # first deriv of kernel
    d2K_array::Array{Float64, 3} = Array{Float64}(undef, 1, 1, 1) # second deriv of kernel
    K_curr::Array{Float64, 3} = Array{Float64}(undef, 1, 1, 1) # kernel values for frequencies in ω_array
    K_new::Array{Float64, 3} = Array{Float64}(undef, 1, 1, 1) # temporary kernel values for frequencies in updated ω_array

    Gbar::Array{Float64, 1} = Array{Float64}(undef, 1) # G(τ) for SAC spectrum
    Gbar_new::Array{Float64, 1} = Array{Float64}(undef, 1) # temporary G(τ) for updated SAC spectrum
    χ2::Float64 = 0. # goodness-of-fit for SAC spectrum
    χ2_new::Float64 = 1.0e20 # temporary goodness-of-fit for updated SAC spectrum
    χ2_min::Float64 = 1.0e20 # minimum goodness-of-fit for configuration within an annealing run
    
    accept_rates::Array{Float64, 3} = Array{Float64}(undef, 1, 1, 1) # acceptance rates of various frequency updates
    amplitude_accept_rates::Array{Float64, 2} = Array{Float64}(undef, 1, 1) # acceptance rates of various amplitude updates
    

    ∆ω_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1) # frequendy update windows
    n_multiδ::Array{Int64, 2} = Array{Int64}(undef, 1, 1) # number of deltas useed in the multi delta edge updates
    
    sampled_spec::Array{Float64, 3} = Array{Float64}(undef, 1, 1, 1) # accumalted spectrum 
    

    χ2_res::Array{Float64, 1} = zeros(Float64, 3) # accumlated χ2
    edge_res::Array{Float64, 2} = zeros(Float64, (3, 2)) # accumlated edge loc.
    n0_res::Array{Float64, 2} = zeros(Float64, (3, 2)) # accumlated n0 value
    ωn0_res::Array{Float64, 2} = zeros(Float64, (3, 2)) # accumlated frquency location of n0

    ω_avg::Array{Float64, 2} = Array{Float64}(undef, 1, 1) # average values of ω_array
    A_avg::Array{Float64, 2} = Array{Float64}(undef, 1, 1) # average values of A_array

    # Temporary arrays used in various updating/measurement steps
    temp_ω_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1) 
    temp_∆ω_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1)
    temp_n_multiδ::Array{Int64, 2} = Array{Int64}(undef, 1, 1)
    temp_A_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1)
    temp_d_array::Array{Float64, 1} = Array{Float64}(undef, 1)
    multiδ_ω_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1)

    # Saved configuration used to reset annealing run after identifying χ2_min
    saved_ω_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1)
    saved_∆ω_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1)
    saved_n_multiδ::Array{Int64, 2} = Array{Int64}(undef, 1, 1)
    saved_A_array::Array{Float64, 2} = Array{Float64}(undef, 1, 1)


    # Output filenames
    anneal_file::String
    accept_rate_file_1::String # rightward decaying edge
    accept_rate_file_2::String # leftward decaying edge
    sample_file::String
    error_file::String
    output_folder::String

    pass::Int64 = 0 # variable used in edge frequency update

    R_L_index::Array{Int64, 1} = [] # [1] if  single_edge or double_edge_symm, [1, 2] else

    χ2_anneal::Array{Float64, 1} = Array{Float64}(undef, 1) # list of χ2 during annealing runs

    fix_edge::Int64 = 0 # if 0 edge is sampled,
                      # if 1 edge is fixed to ω_floor (only for single edge case)
    kernel_type::Symbol # finiteT or zeroT or bosonic
    
    mode::Symbol
    
end
end
SAC = SACs.SAC

###########################################################################################
## Intialization Functions
###########################################################################################

# Read in info from G_file
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
    
    ω_window = log(1/G[N_τ])/τ[N_τ]
    
    self.edge_guess = -(log(G[end])-log(G[end-1]))/(τ[end] - τ[end-1])
    G_D = transpose(cov) * G


    self.τ = τ
    self.G = G_D
    self.σ_inv = σ_inv
    self.cov = cov
    self.Gbar = zeros(Float64, N_τ)
    self.Gbar_new = zeros(Float64, N_τ)
    return
end

# Initialize kernel
function init_kernal!(self::SAC)
    
    ωi_0 = Int64(floor(self.ω_0/self.δω))
    ωi_m = Int64(ceil(self.ω_m/self.δω))

    N_ωi = Int64(ωi_m-ωi_0+1)
   
    
    K_array = zeros(Float64, (self.N_τ, N_ωi, 2))
    dK_array = zeros(Float64, (self.N_τ, N_ωi, 2))
    d2K_array = zeros(Float64, (self.N_τ, N_ωi, 2))

    K_curr = zeros(Float64, (self.N_τ, self.N_ω+1, 2))
    K_new = zeros(Float64, (self.N_τ, self.N_ω+1, 2))

    K_vals = zeros(Float64, self.N_τ)

    if self.kernel_type == :finiteT
        K_function = finiteT_K
    elseif self.kernel_type == :zeroT
        K_function = zeroT_K
    elseif self.kernel_type == :bosonic
        K_function = bosonic_K
    else
        throw("Invalid Kernel type.")
    end

    for i=ωi_0:ωi_m
        ω = i * self.δω
        

        if self.mode == :single_edge 
           
            K_vals .= K_function.(ω, self.τ, self.β) 
            K_array[:, i+1-ωi_0, 1] .= K_vals
            K_array[:, i+1-ωi_0, 1] .= transpose(self.cov) * K_array[:, i+1-ωi_0, 1]
        
        elseif self.mode == :double_edge_in
            K_vals .= K_function.(ω, self.τ, self.β)
            K_array[:, i+1-ωi_0, 1] .= K_vals
            K_array[:, i+1-ωi_0, 1] .= transpose(self.cov) * K_array[:, i+1-ωi_0, 1]

            K_vals .= K_function.(-ω, self.τ, self.β)
            K_array[:, i+1-ωi_0, 2] .= K_vals
            K_array[:, i+1-ωi_0, 2] .= transpose(self.cov) * K_array[:, i+1-ωi_0, 2]
        
        elseif self.mode == :double_edge_out
            K_vals .= K_function.(ω, self.τ, self.β)
            K_array[:, i+1-ωi_0, 1] .= K_vals
            K_array[:, i+1-ωi_0, 1] .= transpose(self.cov) * K_array[:, i+1-ωi_0, 1]

            K_vals .=  K_function.(-ω, self.τ, self.β)
            K_array[:, i+1-ωi_0, 2] .= K_vals
            K_array[:, i+1-ωi_0, 2] .= transpose(self.cov) * K_array[:, i+1-ωi_0, 2]
        
        elseif self.mode == :double_edge_symm
            K_vals .= K_function.(ω, self.τ, self.β) + K_function.(-ω, self.τ, self.β)
            K_array[:, i+1-ωi_0, 1] .= K_vals
            K_array[:, i+1-ωi_0, 1] .= transpose(self.cov) * K_array[:, i+1-ωi_0, 1]

        end
    end
    for R_L in self.R_L_index
        for i=ωi_0+1:ωi_m - 1
            dK_array[:, i+1-ωi_0, R_L] .=  (K_array[:, i+2-ωi_0, R_L] - K_array[:, i-ωi_0, R_L]) / (2 * self.δω)
            d2K_array[:, i+1-ωi_0, R_L] .=  (K_array[:, i+2-ωi_0, R_L] - (2 * K_array[:, i+1-ωi_0, R_L]) + K_array[:, i-ωi_0, R_L]) / (2 * self.δω^2)
        end
        dK_array[:, 1, R_L] =  (K_array[:, 2, R_L] - K_array[:, 1, R_L]) / self.δω
        dK_array[:, ωi_m+1, R_L] = (K_array[:, ωi_m+1, R_L] - K_array[:, ωi_m, R_L]) / self.δω
    end

    self.K_array, self.dK_array, self.d2K_array = K_array,  dK_array, d2K_array
    self.K_curr, self.K_new = K_curr, K_new
    self.ωi_0, self.ωi_m, self.N_ωi = ωi_0, ωi_m, N_ωi

end

# Funtion to calculate kernel that can better handle large β
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


# Initialize arrays in struct and outfiles
function init_struct!(self::SAC)


    self.ω_array = zeros(Float64, (self.N_ω + 1, 2))
    self.d_array = zeros(Float64, (self.N_ω + 1, 2))
    self.A_array = zeros(Float64, (self.N_ω + 1, 2))

    self.n_multiδ = zeros(Int64, (self.N_ω + 1, 2))
    self.∆ω_array = zeros(Float64, (self.N_ω + 1, 2))
    self.accept_rates = zeros(Float64, (2, self.N_ω + 1, 2))
    self.amplitude_accept_rates = zeros(Float64, (3, 2))


    self.sampled_spec = zeros(Float64, (self.N_ωi, 2, 2))
    
    self.ω_avg = zeros(Float64, (self.N_ω + 1, 2))
    self.A_avg = zeros(Float64, (self.N_ω + 1, 2))

    self.χ2_res = zeros(Float64, 3)
    self.edge_res = zeros(Float64, (3, 2))
    self.n0_res = zeros(Float64, (3, 2))
    self.ωn0_res = zeros(Float64, (3, 2))

    self.temp_ω_array = copy(self.ω_array) 
    self.temp_∆ω_array = copy(self.∆ω_array)
    self.temp_n_multiδ = copy(self.n_multiδ)
    self.temp_A_array = copy(self.A_array)
    self.temp_d_array = zeros(Float64, self.N_ω + 1)

    self.multiδ_ω_array = copy(self.ω_array)

  

    self.χ2_anneal = zeros(Float64, self.N_anneal)


    # column names should be self-evident
    # wc0,w_cf are boundaries of continuum if A_c > 0
    f = open(self.anneal_file, "w")
       println(f, "j,theta,chi2_min,chi2_avg,chi2_sigma,edge_R,wn0_R,n0_R,eps_R,edge_R_sigma,wn0_R_sigma,n0_R_sigma,"*
                                            "edge_L,wn0_L,n0_L,eps_L,edge_L_sigma,wn0_L_sigma,n0_L_sigma,"*
                                            "wc0,w_cf")
    close(f)

    #acceptance rates:
        # continuum freq. and update window
        # continuum amp.
        # single edge freq. and update window
        # multi edge freq. and avearge cluster size
        # edge amp
    f = open(self.accept_rate_file_1, "w")
       println(f, "j,ar_cont,∆w_cont,ar_cont_amp,ar_edge_s,∆w_edge,ar_edge_m,n_clust,ar_edge_amp")
    close(f)

    if self.mode != :single_edge || self.mode != :double_edge_symm
        f = open(self.accept_rate_file_2, "w")
           println(f, "j,ar_cont,∆w_cont,ar_cont_amp,ar_edge_s,∆w_edge,ar_edge_m,n_clust,ar_edge_amp")
        close(f)
    end


    f = open(self.sample_file, "w")
       println(f, "j,theta,chi2_min,chi2_avg,chi2_sigma,edge_R,wn0_R,n0_R,eps_R,edge_R_sigma,wn0_R_sigma,n0_R_sigma,"*
                                            "edge_L,wn0_L,n0_L,eps_L,edge_L_sigma,wn0_L_sigma,n0_L_sigma,"*
                                            "wc0,w_cf")
    close(f)

end

# Setting amplitudes in A_array
# Refer to https://www.sciencedirect.com/science/article/pii/S0370157322003921 for form of A_e[i]
function set_amplitudes!(self::SAC, m::Int64, R_L::Int64)
    if m == 0
        self.n0 = [self.N_e, self.N_e] .* 0.75
        self.∆n0 .= 2.
        

        self.ε_0 = [0.5, 0.5]
        self.∆ε_0 .= 0.05
    end

    if m == -1
        self.∆n0 .= 2.
        self.∆ε_0 .= 0.1
    end
    self.A_array[1, R_L] = 0
    self.temp_A_array[1, R_L] = 0

    ε_1 = (log(2) * self.c[R_L])^2 * self.ε_0[R_L]
    for i=1:self.N_e
        x = log(i/self.n0[R_L])
        if self.c[R_L] < 0
            ln_Ai = self.c[R_L]*x + sqrt((self.c[R_L]*x)^2 + ε_1)
        else
            ln_Ai = self.c[R_L]*x - sqrt((self.c[R_L]*x)^2 + ε_1)
        end

        self.A_array[i+1, R_L] = exp(ln_Ai / 2)
    end

    # Continuum deltas with equal amplitudes
    if self.N_ω > self.N_e
        for i=(self.N_e+1):self.N_ω
            self.A_array[i+1, R_L] = 1.
        end
    end    

    # Only rightward delta ω array contains continuum
    f1, f2 = 0., 0.
    if (R_L == 1) 
        f1 = self.A_r / (1 + self.A_c) # weight of rightward decaying edge
        f2 = self.A_c / (1 + self.A_c) # weight of continuum 
    else
        f1 = (1-self.A_r) / (1 + self.A_c) # weight of leftward decaying edge
        f2 = 0.
    end

    self.A_array[2:self.N_e+1, R_L] .*= (f1 / sum(self.A_array[2:self.N_e+1, R_L]))
    
    self.A_array[(self.N_e+2):end, R_L] .*= (f2 / sum(self.A_array[(self.N_e+2):end, R_L]))
    
    self.temp_A_array[:, R_L] .= self.A_array[:, R_L]

    
end


# Setting initial configuration of delta functions
# This step is crucial, as it is challenging for the program to equillibrate 
# if the intial config is far from optimal
function init_config_jk!(self::SAC, j::Int64, k::Int64)



    ∆d = 0.


    """ 
    If you run into equillibration issues, i.e. can't find good initial config so χ2 is
    remains very large during anneal, then change parameters used to set self.ω_array[2, :]
    and ∆d. They were chosen by trial and error and may depend on specifics of 
    spectrum that you are trying to resolve.

    Specifically the denominators in (j*(self.ω_m)/~denom~) and ((k+5)/~denom~)
    """


    k0 = 55

    if self.fix_edge == 1
        self.ω_array[2, :] .= self.ω_array[1, :] 
        ∆d = abs(((self.ω_m * (k+5)/k0) - self.ω_array[2, 1]))

    elseif self.mode == :single_edge
        self.ω_array[2, 1] = self.ω_array[1, 1] + (sign(self.edge_guess) * ((j-1)*(self.ω_m)/500))
        ∆d = abs(((self.ω_m * (k+5)/k0) - self.ω_array[2, 1]))

    elseif self.mode == :double_edge_in
        self.ω_array[2, 1] = self.ω_array[1, 1] + (sign(self.edge_guess) * (j*(self.ω_m)/500))
        self.ω_array[2, 2] = self.ω_array[2, 1] - ((k+5)*(self.ω_m)/k0)

        # seting distances between deltas so fit between two edges
        ∆d = abs((-self.ω_array[2, 2] - self.ω_array[2, 1])) * 0.5
    
    elseif self.mode == :double_edge_out
        self.ω_array[2, 1] = self.ω_array[1, 1] + (j*(self.ω_m)/500)
        self.ω_array[2, 2] = self.ω_array[1, 2] + (k*(self.ω_m)/500)

        # ∆d = abs(((self.ω_m * (k+5)/k0) - self.ω_array[2, 1]))
        ∆d = abs(((self.ω_m * (0.5)) - self.ω_array[2, 1]))

    elseif self.mode == :double_edge_symm
        self.ω_array[2, 1] = self.ω_array[1, 1] + (j*(self.ω_m)/500)
        ∆d = abs(((self.ω_m * (k+5)/k0) - self.ω_array[2, 1]))
    end


    self.d_array[2:self.N_e+1, :] .= sqrt.(1:self.N_e)
    self.d_array[3:self.N_e+1, :] .*= ∆d/sum(self.d_array[3:self.N_e+1, :])
    
    # edge deltas
    for i=2:self.N_e
        for R_L in self.R_L_index
            self.ω_array[i+1, R_L] = self.ω_array[i, R_L] + self.d_array[i+1, R_L]
        end
    end

    # continuum deltas
    if self.N_e < self.N_ω
        self.ω_array[self.N_e+2:end, 2] .= self.ω_array[3, 2]

        self.ω_array[self.N_e+2, 1] = self.ω_array[3, 1]
        if self.mode == :double_edge_in
            ∆d = abs((-self.ω_array[2, 2] - self.ω_array[self.N_e+2, 1])) * 0.5
        else
            ∆d = abs((self.ω_array[self.N_e+1, 1] - self.ω_array[self.N_e+2, 1])) * 0.5
        end

        self.d_array[self.N_e + 2:end] .= ∆d/(self.N_ω - self.N_e)
        
        for i=self.N_e+2:self.N_ω
            self.ω_array[i+1, 1] = self.ω_array[i, 1] + ∆d/(self.N_ω - self.N_e)
        end
    end


    # println([self.edge_guess, self.ω_array[2, 1], maximum(self.ω_array[2, 1])])
    # for i=2:self.N_e+1
    #     println(self.ω_array[i, 1])
    # end
  
end



# Loop over spectra with different edge location and different `widths', 
# setting the initial spectrum to the best one found.
function init_config_dual!(self::SAC)

    for R_L in self.R_L_index
        set_amplitudes!(self, 0, R_L)
    end

    j_max = 50
    k_max = 50

    
    # if edge is fixed, set ω_array[1, :] to ω_floor, otherwise set to 0
    # and the sign of edge_guess (from G(τ)) will tell whether to scan in + or - freq. dir.

    if self.fix_edge == 1
        self.ω_array[1, :] .= self.ω_floor
        j_max = 1
    else
        self.ω_array[1, :] .= 0.
    # elseif self.mode == :double_edge_in       
    #     self.ω_array[1, :] .= 0.
    # elseif self.mode == :double_edge_out   
    #     self.ω_array[1, :] .= 0#self.edge_guess
    # else
    #     self.ω_array[1, :] .= 0. #self.edge_guess
    #     # self.ω_array[:, 2] .= 0.
    #     # j_max = 1
    end
    
    self.χ2_min = 1e20  

  


    k_opt, j_opt = 1, 1
    for k=1:k_max
        for j=1:j_max
            init_config_jk!(self, j, k)
            
            if self.mode == :double_edge_in
                # Left tail error, skip (j, k) pair
                if self.ω_array[2, 1] > - self.ω_array[self.N_e+1, 2]
                    continue
                end

                # Right tail error, skip (j, k) pair
                if maximum(self.ω_array[:, 1]) > - self.ω_array[2, 2]
                    continue
                end
            end

        
            calc_Gbar!(self, 1)
            
            χ2 = calc_χ2(self)
            # println([χ2, self.ω_array[2, 1], maximum(self.ω_array[:, 1]), self.ω_m * (k+5)/105])
            if χ2 < self.χ2_min
                k_opt = k
                j_opt = j
                self.χ2_min = χ2
            end
        end
    end


    init_config_jk!(self, j_opt, k_opt)

    write_log(self, "Rightward Edge: [$(self.ω_array[2, 1]), $(self.ω_array[self.N_e+1, 1])] ($j_opt, $k_opt)")
    
    if self.mode == :double_edge_in || self.mode == :double_edge_out
        write_log(self, "Leftward Edge: [$(-self.ω_array[2, 2]), $(-self.ω_array[self.N_e+1, 2])]")
    end
    
    calc_Gbar!(self, 1)
    
    self.χ2 = calc_χ2(self)

    if self.A_c == 0
        @assert self.χ2 ≈ self.χ2_min 
    end

    for R_L in self.R_L_index
        self.n_multiδ[2:self.N_e - 1, R_L] .= 1 + (self.N_ω ÷ 20)
        for i=2:self.N_e-1
            self.n_multiδ[i, R_L] = min(self.n_multiδ[i, R_L], self.N_e - i - 1)
        end
        self.n_multiδ[self.N_e + 2:end, R_L] .= -1


        self.∆ω_array[2, R_L] = self.d_array[3, R_L]
        for i=2:self.N_ω
            self.∆ω_array[i+1, R_L] = self.d_array[i+1, R_L]
        end
        check_ω(self, 0., R_L)
    end


end


###########################################################################################
## Utility Functions
###########################################################################################

# Calculate G(τ) for delta config
function calc_Gbar!(self::SAC, update::Int64)
    if update == 1
        for R_L in self.R_L_index
            @inbounds @simd for i = 1:self.N_τ
                ω = floor(Int64, self.ω_array[1, R_L]/self.δω)
                d = self.ω_array[1, R_L] - (ω * self.δω)
                self.K_curr[i, 1, R_L] = (self.K_array[i, ω+1-self.ωi_0, R_L] +
                                    (self.dK_array[i, ω+1-self.ωi_0, R_L]*d) +
                                    (self.d2K_array[i, ω+1-self.ωi_0, R_L]*d^2))

                @inbounds @simd for j=1:self.N_ω
                    ω = floor(Int64, self.ω_array[j+1, R_L]/self.δω)
                    d = self.ω_array[j+1, R_L] - (ω * self.δω)
                    self.K_curr[i, j+1, R_L] = (self.K_array[i, ω+1-self.ωi_0, R_L] +
                                    (self.dK_array[i, ω+1-self.ωi_0, R_L]*d) +
                                    (self.d2K_array[i, ω+1-self.ωi_0, R_L]*d^2))
                end
            end
        end
        

        
    end
    # add rightward edge contribution
    self.Gbar .= (self.K_curr[:, :, 1] * self.A_array[:, 1])

    # add leftward edge contribution unless single edge or symmetric double edge
    if self.mode == :double_edge_in || self.mode == :double_edge_out 
        self.Gbar .+= (self.K_curr[:, :, 2] * self.A_array[:, 2])
    end

end

# Calculation G(τ) after changing kth delta to ω_k
function calc_Gbar_new!(self::SAC, k::Int64, ω_k::Float64, R_L::Int64)
    ω = floor(Int64, ω_k/self.δω)
    d = (ω_k - (ω * self.δω))
    d2 = d^2
    A_k = self.A_array[k, R_L]
    @inbounds @simd for i=1:self.N_τ
        self.K_new[i, k, R_L] = (self.K_array[i, ω+1-self.ωi_0, R_L] + 
                           (self.dK_array[i, ω+1-self.ωi_0, R_L]*d) + 
                           (self.d2K_array[i, ω+1-self.ωi_0, R_L]*d2))
        self.Gbar_new[i] += A_k * (self.K_new[i, k, R_L] - self.K_curr[i, k, R_L])
    end
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


# Run series of tests to ensure deltas spaced properly and in correct boundaries
function check_ω(self::SAC, m::Float64, R_L::Int64)

    #intialize dummy variables to use for various ω_array tests
    var1 = 0.
    var2 = 0.

    if self.ω_array[2, R_L] < self.ω_floor[R_L]
        f = open(self.error_file,  "w")
            println(f, m, " ", R_L,  " Leading delta violation: ", self.ω_floor[R_L], " ", self.ω_array[1, R_L], " ", self.ω_array[2, R_L])
        close(f)
        exit() 
    end


    if minimum(self.ω_array[:, R_L]) < self.ω_floor[R_L] || maximum(self.ω_array[:, R_L]) > self.ω_m
        f = open(self.error_file,  "w")
            println(f, m, " ", R_L,  " Out-of-bound frequency: ", i+1, " ", minimum(self.ω_array[:, R_L]), " ", maximum(self.ω_array[:, R_L]) ," ",self.ω_floor[R_L], " ", self.ω_m)
        close(f)
        exit() 
    end
 


    if self.N_e<self.N_ω && R_L == 1
        var1 = minimum(self.ω_array[self.N_e + 2:end, R_L])
        if self.ω_array[2, R_L] > var1
            f = open(self.error_file,  "w")
                println(f, m, " ", R_L,  " Lower continuum bound violation: ", var1, " ", self.ω_array[2, R_L])
            close(f)
            exit() 
        end
    end

    var1 = self.ω_array[3, R_L] - self.ω_array[2, R_L]
    if var1 < self.d0
        f = open(self.error_file,  "w")
            println(f, m, " ", R_L,  " Min separation violation: ", var1, " ", self.d0)
        close(f)
        exit() 
    end

    for i=2:self.N_e-1
        var1 = self.ω_array[i+1, R_L] -  self.ω_array[i, R_L]
        var2 = self.ω_array[i+2, R_L] -  self.ω_array[i+1, R_L]
        if var2 < var1 
            if var1 - var2 > 1e-10
                f = open(self.error_file,  "w")
                    println(f, m, " ", R_L,  " Order violation: ", "delta_number", " ","delta_location", " ", "d_minus", " ", "d_plus")
                    println(f, m, " ", R_L,  " Order violation: ", i+1, " ",self.ω_array[i+1, R_L], " ", var1, " ", var2)
                close(f)
                exit() 
            end

        end
    end

    f = open(self.error_file,  "w")
        println(f, m, " ", R_L,  " No Errors")
    close(f)

end

# Check if correct spacing is maintained after changing kth delta to ω_k
function check_order(self::SAC, k::Int64, ω_k::Float64, R_L::Int64)
    d1, d2, d3, d4 = -1, -1, -1, -1
    if k == 2

        d1 = self.ω_array[3, R_L] -  ω_k
        if d1 < self.d0
            return 0
        end
        if self.N_e > 2
            d2 = self.ω_array[4, R_L] - self.ω_array[3, R_L]
            if d2 < d1 
                return 0
            end
        end
    elseif k==3
        d1 = ω_k - self.ω_array[2, R_L]
        if d1<self.d0
            return 0
        end
        if self.N_e > 2
            d2 = self.ω_array[4, R_L] - ω_k
            if d2 < d1 
                return 0
            end
        end
        if self.N_e > 3
            d3 = self.ω_array[5, R_L] - self.ω_array[4, R_L]
            if d3 < d2
                return 0
            end
        end
    elseif k == self.N_e
        d1 = self.ω_array[k-1, R_L] - self.ω_array[k-2, R_L]
        d2 = ω_k - self.ω_array[k-1, R_L] 
        if d2 < d1 
            return 0
        end
        d3 = self.ω_array[k+1, R_L] - ω_k
        if d3 < d2
            return 0
        end
    elseif k == self.N_e+1
        d1 = self.ω_array[k-1, R_L] - self.ω_array[k-2, R_L]
        d2 = ω_k - self.ω_array[k-1, R_L]
        if d2 < d1 
            return 0
        end
    else
        d1 = self.ω_array[k-1, R_L] - self.ω_array[k-2, R_L]
        d2 = ω_k - self.ω_array[k-1, R_L]
        if d2 < d1 
            return 0
        end
        d3 = self.ω_array[k+1, R_L] - ω_k
        if d3 < d2
            return 0
        end
        d4 =  self.ω_array[k+2, R_L] - self.ω_array[k+1, R_L]
        if d4 < d3
            return 0
        end
    end
    return 1
end



##########################################################################################
# Update Functions
##########################################################################################

# Update freq. of single delta in edge
function single_edge_move(self, θ::Float64, R_L)
 

    if self.A_c == 1
        return nothing
    end

    if self.mode == :double_edge_in
        # negative sign because the leftward deltas are stored as their negatives
        self.ω_b0[R_L] = - maximum(self.ω_array[2:end, (3-R_L)]) #first delta cant go to right(R_L=1)/left(R_L=2) of this value (tail)
        self.ω_bf[R_L] = - self.ω_array[2, (3-R_L)] #last delta cant go to right(R_L=1)/left(R_L=2) of this value (edge)
    else
        # no restrictions if single_edge or double_edge_out
        self.ω_b0[R_L] = self.ω_m
        self.ω_bf[R_L] = self.ω_m
    end

    # if there is a continuum, rightward delta edge cant go past lowest delta in continuum
    if self.A_c > 0 && R_L == 1
        self.ω_b0[R_L] = min(minimum(self.ω_array[self.N_e + 2:end, 1]), self.ω_b0[R_L])
    end



    k::Int64 = 0
    ω_k::Float64 =0.
    χ2_new::Float64 = 0.
    P::Float64 = 0.

    k_0::Int64 = 2

    # if edge is fixed, don't update first delta
    if self.fix_edge == 1
        k_0 = 3
    end

    for i=1:self.N_e

        if mod(i, 4) == 0 && self.fix_edge == 0
            k = 2
        else
            k = rand(k_0:self.N_e+1)
        end
    

        ω_k = self.ω_array[k, R_L] + (self.∆ω_array[k, R_L] * (rand()-0.5))

        if k == 2
            if ω_k > self.ω_b0[R_L] 
                continue
            end
        end

        if k < self.N_e+1
            if ω_k < self.ω_array[k-1, R_L] || ω_k > self.ω_array[k+1, R_L]
                continue
            end
        else
            if ω_k < self.ω_array[k-1, R_L] || ω_k > self.ω_m || ω_k > self.ω_bf[R_L] #tail cant go past edge
                continue
            end 
        end

     

        self.pass = check_order(self, k, ω_k, R_L)
        if self.pass == 1
            self.Gbar_new .= self.Gbar
            calc_Gbar_new!(self, k, ω_k, R_L)
            χ2_new = calc_χ2_new(self)
            P = exp((self.χ2 - χ2_new)/(2*θ))
            if rand() < P

        
                @inbounds @simd for j=1:self.N_τ
                    self.Gbar[j] = self.Gbar_new[j]
                    self.K_curr[j, k, R_L] = self.K_new[j, k, R_L]
                end

                
                self.ω_array[k, R_L] = ω_k
                self.χ2 = χ2_new
                if χ2_new < self.χ2_min
                    self.χ2_min = χ2_new
                end
               self.accept_rates[1, k, R_L] += 1
            end
        else
        end
    end
    return nothing

end


# Update freq. of multiple deltas in edge
function multi_edge_move(self::SAC, θ::Float64, R_L::Int64)
    
    if self.A_c == 1
        return nothing
    end

    n::Int64 = 0
    k::Int64 = 0
    j::Int64 = 0.
    ω_k::Float64 =0.
    ω_j::Float64 =0.
    χ2_new::Float64 = 0.
    P::Float64 = 0.

    if self.N_e < 4
        return nothing
    end

    if self.mode == :double_edge_in
        #negative because the leftward deltas are stored as their negatives
        self.ω_b0[R_L] = - maximum(self.ω_array[:, (3-R_L)]) #first delta cant go to right(R_L=1)/left(R_L=2) of this value
        self.ω_bf[R_L] = - self.ω_array[2, (3-R_L)] #last delta cant go to right(R_L=1)/left(R_L=2) of this value
    else
        self.ω_b0[R_L] = self.ω_m
        self.ω_bf[R_L] = self.ω_m
    end

    if self.A_c > 0 && R_L == 1
        self.ω_b0[R_L] = min(minimum(self.ω_array[self.N_e + 2:end, 1]), self.ω_b0[R_L])
    end

   
    for i=0:((self.N_e -1))
        # update edge location more frequently, 
        # uses slightly different multi delta update method
        if mod(i, 4) == 0 && self.fix_edge == 0
            k = 2
            n = self.n_multiδ[2, R_L]
            multi_1!(self, n, R_L)

        else
            k = rand(3:self.N_e-1)
            n = self.n_multiδ[k, R_L]
            multi_k!(self, k, n, R_L)
        end

        if self.pass == 1
            self.Gbar_new .= self.Gbar
            ω_js = self.multiδ_ω_array[k:k+n-1, R_L]
            for j=1:n
                calc_Gbar_new!(self, k + j - 1, ω_js[j] , R_L)
            end

            χ2_new = calc_χ2_new(self)
            P = exp((self.χ2 - χ2_new)/(2*θ))
            if rand() < P
          
                @inbounds @simd for t=1:self.N_τ
                    self.Gbar[t] = self.Gbar_new[t]
                     @inbounds @simd for j=k:k+n-1
                        self.K_curr[t, j, R_L] = self.K_new[t, j, R_L]
                    end
                end

                @inbounds @simd for j=k:k+n-1
                    self.ω_array[j, R_L] = self.multiδ_ω_array[j, R_L]
                end
                
                self.χ2 = χ2_new
                if χ2_new < self.χ2_min
                    self.χ2_min = χ2_new
                end
                
                self.accept_rates[2, k, R_L] += 1
                
            end
        end
    end
    return nothing
end

# Fuction used for multi delta update of edge 
# moves deltas around while ensuring that max distance doesn't exceed
# that for delta right above cluster
function multi_1!(self::SAC, n::Int64, R_L::Int64)
    self.pass = 0

    d1 = self.ω_array[n+2, R_L] - self.ω_array[n+1, R_L]
    d_array = view(self.temp_d_array, 1:n-1)
    @inbounds @simd for i=2:n
        d_array[i-1] = self.d0 + ((d1 - self.d0) * rand())
    end

    ω_1 = self.ω_array[n+1, R_L] - sum(d_array)
  


    if ω_1 > self.ω_floor[R_L] && ω_1 < self.ω_b0[R_L] #edge must be below tail
        if n > 2
            sort!(d_array)
        end
        self.multiδ_ω_array[2, R_L] = ω_1
        @inbounds @simd for i=2:n
            self.multiδ_ω_array[i+1, R_L] = self.multiδ_ω_array[i, R_L] + d_array[i-1]
        end
    
        self.pass = 1
    end

    return nothing
end

# Fuction used for multi delta update for deltas above edge 
# moves deltas around while ensuring that max distance doesn't exceed
# that for delta right above cluster, using more sophisticated method
# refer to https://www.sciencedirect.com/science/article/pii/S0370157322003921 for details
function multi_k!(self::SAC, k::Int64, n::Int64, R_L::Int64)
    self.pass = 0

    i::Int64 = 0
    j::Int64 = 0
    di::Float64 = 0.
    dj::Float64 = 0.
    
    
    d_min::Float64 = 0.
    d_max::Float64 = 0.
    d_array = view(self.temp_d_array, 1:n)
   
    @inbounds @simd for i=0:n-1
       d_array[i+1] = self.ω_array[k+i+1, R_L] - self.ω_array[k+i, R_L]
    end

    
    dk_minus1::Float64 =  self.ω_array[k, R_L] - self.ω_array[k-1, R_L]
    dk_plusn::Float64 = self.ω_array[k+n+1, R_L] - self.ω_array[k+n, R_L]
    
    for _=1:(n ÷ 2)
        i = rand(0:n-1)
        j = i
        while j == i
            j = rand(0:n-1)
        end
        
        di = d_array[i+1]
        dj = d_array[j+1]
        
        d_min = max(dk_minus1, di + dj - dk_plusn)
        d_max = min(dk_plusn, di + dj - dk_minus1)
        
        d_array[i+1] = d_min + (d_max - d_min) * rand()
        d_array[j+1] = di + dj - d_array[i+1]
    end
    sort!(d_array)
    self.multiδ_ω_array[k, R_L] = self.ω_array[k, R_L]
    @inbounds @simd for i=k+1:k+n-1
        self.multiδ_ω_array[i, R_L] = self.multiδ_ω_array[i-1, R_L] + d_array[i - k]
    end

    if self.multiδ_ω_array[k+n-1, R_L] < self.ω_bf[R_L]
        self.pass = 1
    end
    return nothing
end


# Function for updating amplitudes of edge deltas
# just changes the value of n0 and ε, so will have no effect if p = 1/2
function ampltiude_edge_move(self::SAC, θ::Float64, R_L::Int64)
    
    if self.A_c == 1
        return nothing
    end

    n1::Float64 = self.n0[R_L]
    ε_1::Float64 = self.ε_0[R_L]
    ln_Ai::Float64 = 0

    if self.N_e < 3
        return
    end

    ε = (log(2) * self.c[R_L])^2 
   
    n1 += (self.∆n0[R_L] * (rand() - 0.5))
    if n1 < 5 || n1 > (self.N_e - 5)
        return
    end

    # No rounding, unclear if this has a major effect 
    # ε_0_prime = 0
    # self.∆ε_0[R_L] = 0

    # Rounding 
    ε_1 += self.∆ε_0[R_L] * (rand() - 0.5)

    if ε_1 < 0 || ε_1 > 1
        return
    end


    
    for i=1:self.N_e
        x = log(i/n1)
        if self.c[R_L] < 0
            ln_Ai = self.c[R_L]*x + sqrt((self.c[R_L]*x)^2 + ε_1*ε)
        else
            ln_Ai = self.c[R_L]*x - sqrt((self.c[R_L]*x)^2 + ε_1*ε)
        end
        self.temp_A_array[i+1, R_L] = exp(ln_Ai/2)

    end
    
    f1 = 0.
    if (R_L == 1) 
        f1 = self.A_r / (1 + self.A_c)
    else
        f1 = (1-self.A_r) / (1 + self.A_c)
    end


    self.temp_A_array[2:self.N_e+1, R_L] *= (f1 / sum(self.temp_A_array[2:self.N_e+1, R_L]))
    

    self.Gbar_new .= (self.Gbar .+ 
                      self.K_curr[:, 2:self.N_e+1, R_L] *
                      (self.temp_A_array[2:self.N_e+1, R_L] .- self.A_array[2:self.N_e+1, R_L]))



    χ2_new = calc_χ2_new(self)
    P = exp((self.χ2 - χ2_new)/(2*θ))
    if rand() < P
        self.n0[R_L] = n1
        self.ε_0[R_L] = ε_1

        self.A_array[2:self.N_e+1, R_L] .= self.temp_A_array[2:self.N_e+1, R_L]

        @inbounds @simd for j=1:self.N_τ
            self.Gbar[j] = self.Gbar_new[j]
            
        end

        self.χ2 = χ2_new
        if χ2_new < self.χ2_min
            self.χ2_min = χ2_new
        end


        self.amplitude_accept_rates[1, R_L] += 1

    end
   return nothing
end

# Function for updating freq. of continuum deltas
# for R_L = 1 only
function single_cont_move(self::SAC, θ::Float64) 
    
    if self.N_e == self.N_ω || self.A_c == 0
        return nothing
    end

    # continuum must lie within bounds of edge spectrum
    ωc_min = self.ω_array[2, 1] 
    ωc_max = self.ω_array[self.N_e + 1, 1]
    
    k::Int64 = 0
    ω_k::Float64 =0.
    χ2_new::Float64 = 0.
    P::Float64 = 0.

    for i=(self.N_e+2):self.N_ω + 1

        k = rand((self.N_e+2):self.N_ω+1)
        ω_k = self.ω_array[k] + (self.∆ω_array[k, 1] * (rand()-0.5))

        if ω_k < ωc_min || ω_k > ωc_max
            continue
        end

        self.Gbar_new .= self.Gbar
        calc_Gbar_new!(self, k, ω_k, 1)
        χ2_new = calc_χ2_new(self)
        P = exp((self.χ2 - χ2_new)/(2*θ))
        if rand() < P
            @inbounds @simd for j=1:self.N_τ
                self.Gbar[j] = self.Gbar_new[j]
                self.K_curr[j, k, 1] = self.K_new[j, k, 1]
            end

            
            self.ω_array[k, 1] = ω_k
            self.χ2 = χ2_new
            if χ2_new < self.χ2_min
                self.χ2_min = χ2_new
            end
           self.accept_rates[1, k, 1] += 1
        end 
    end
end

# Function for updating amp. of continuum deltas
# conserves first moment of two updated deltas so that A_c doesn't change 
# for R_L = 1 only
function ampltiude_cont_move(self::SAC, θ::Float64)

    if self.N_e == self.N_ω || self.A_c == 0
        return
    end

    N_accept::Int64 = 0
    k::Int64 = 0
    j::Int64 = 0
    A_j::Float64 = 0
    A_k::Float64 = 0
    
    for i=1:self.N_ω
        j = rand((self.N_e+2):self.N_ω)
        k = j
        while k == j
            k = rand((self.N_e+2):self.N_ω)
        end

        A_j = self.A_array[j, 1]
        A_k = self.A_array[k, 1]
        
        m0 = A_j + A_k
        
        r = rand()

        δA_j = (r * m0) - A_j
        δA_k = ((1 - r) * m0) - A_k
        
        

        @inbounds @simd for t=1:self.N_τ
            self.Gbar_new[t] = (self.Gbar[t] + δA_j * self.K_curr[t, j, 1]
                                             + δA_k * self.K_curr[t, k, 1])
        end

        
        χ2_new = calc_χ2_new(self)
        P = exp((self.χ2 - χ2_new)/(2*θ))
        if rand() < P
            self.A_array[j, 1] += δA_j
            self.A_array[k, 1] += δA_k

            @inbounds @simd for t=1:self.N_τ
                self.Gbar[t] = self.Gbar_new[t]
            end

            self.χ2 = χ2_new
            if χ2_new < self.χ2_min
                self.χ2_min = χ2_new
            end

            N_accept += 1
        end
    end

    self.temp_A_array[self.N_e+2:end, 1] .= self.A_array[self.N_e+2:end, 1]
    self.amplitude_accept_rates[2, 1] += N_accept / (self.N_ω - self.N_e)
end


##########################################################################################
## Functions for performing above sampling routines and recording results  
##########################################################################################

# Function to run set of updates at a given θ
function run_updates(self::SAC, steps::Int64, θ::Float64)

    # reset accumalted measureables before running sweep of updates
    self.χ2_res[1] = 0
    self.edge_res[1, :] .= 0
    self.n0_res[1, :] .= 0
    self.ωn0_res[1, :] .= 0
    self.accept_rates .= 0

    
    calc_Gbar!(self, 0)

    for i=1:steps
        for R_L in self.R_L_index

            single_edge_move(self, θ, R_L)
            multi_edge_move(self, θ, R_L)
        
            ampltiude_edge_move(self, θ, R_L)
           
    
        end

         if self.A_c > 0
            single_cont_move(self, θ)

            ampltiude_cont_move(self, θ)
        end
       
        
        measure(self)

    end
    
    # record results for set of sweeps
    record_bin(self, steps)

    for R_L in self.R_L_index
        check_ω(self, θ, R_L)
    end
    

end

# Function to record accumulated spectra and measure various info (edge, n0, etc.)
function measure(self::SAC)
    self.χ2_res[1] += self.χ2
   

    for R_L in self.R_L_index

        self.n0_res[1, R_L] += self.n0[R_L]/self.N_e # fraction of deltas where transition in ampltiude form happens 
        self.ωn0_res[1, R_L] += self.ω_array[floor(Int64, self.n0[R_L]), R_L] # freq. where transition in ampltiude form happens 

        self.A_avg[:, R_L] .+= self.A_array[:, R_L] # accumualte avg amps
        self.ω_avg[:, R_L] .+= self.ω_array[:, R_L] # accumualte avg freqs
       

        ωi_edge = floor(Int64, self.ω_array[2, R_L]/self.δω) # loc. of edge
        self.edge_res[1, R_L] += self.ω_array[2, R_L] # accumualte edfe loc. freq

        # accumulate average spectrum contribution from edge     
        for i=1:self.N_e
            ωi = floor(Int64, self.ω_array[i+1, R_L]/self.δω)
            self.sampled_spec[ωi+1-self.ωi_0, 1, R_L] += self.A_array[i+1, R_L]
        end

        # accumulate average spectrum contribution from continuum  
        for i=(self.N_e+1):self.N_ω
            ωi = floor(Int64, self.ω_array[i+1, R_L]/self.δω)
            self.sampled_spec[ωi+1-self.ωi_0, 2, R_L] += self.A_array[i+1, R_L]
        end
    end
end

# Funtion to accumulate measured info about spectrum in bins, which are used to track 
# fluctuations of the spectrum
function record_bin(self::SAC, steps::Int64)
    self.χ2_res[1] /= steps
    self.n0_res[1, :] ./= steps
    self.ωn0_res[1, :] ./= steps
    self.edge_res[1, :] ./= steps


    self.χ2_res[2] += self.χ2_res[1]
    self.n0_res[2, :] .+= self.n0_res[1, :]
    self.ωn0_res[2, :] .+= self.ωn0_res[1, :]
    self.edge_res[2, :] .+= self.edge_res[1, :]

    self.χ2_res[3] += self.χ2_res[1]^2
    self.n0_res[3, :] .+= self.n0_res[1, :].^2
    self.ωn0_res[3, :] .+= self.ωn0_res[1, :].^2
    self.edge_res[3, :] .+= self.edge_res[1, :].^2
end


# Function to run updating and measurement steps in bins
function run_bins(self::SAC, steps::Int64, bins::Int64, θ::Float64)

    self.χ2_res .= 0
    self.edge_res .= 0
    self.n0_res .= 0
    self.ωn0_res .= 0
    self.sampled_spec .= 0
    self.ω_avg .= 0
    self.A_avg .= 0

  

   for i=1:bins
        run_updates(self, steps, θ)
        for R_L in self.R_L_index
            adjust_∆(self, steps, R_L)
        end

    end
    bin_averages(self, bins, false)


end



# Funciton ot adjust parameters in sampling steps to maintain ~50% acceptance rate
function adjust_∆(self::SAC, steps::Int64, R_L::Int64)

    self.accept_rates[1, :, R_L] ./= steps # single edge freq. update
    self.accept_rates[2, :, R_L] ./= steps # multi edge freq. update


    # adjust freq. window in single delta updates
    for i=1:self.N_ω
        if self.accept_rates[1, i+1, R_L] > 0.55
            self.∆ω_array[i+1, R_L] *= 1.25
        elseif self.accept_rates[1, i+1, R_L] < 0.45
             self.∆ω_array[i+1, R_L] /= 1.2
        end
     end

    
    # adjust cluster size for multi delta updates
    for i=0:self.N_e-2
        if self.accept_rates[2, i+1, R_L] > 0.55
            self.n_multiδ[i+1, R_L] = floor(Int64, self.n_multiδ[i+1, R_L] * 1.25) + 1
        elseif self.accept_rates[2, i+1, R_L] < 0.45
            self.n_multiδ[i+1, R_L] = max(1, floor(Int64, self.n_multiδ[i+1, R_L] / 1.2))
        end
    end

    # ensure none of the cluster sizes are bigger
    # than the number of deltas to the right
    if self.N_e > 3 
        for i=2:self.N_e-1
            self.n_multiδ[i, R_L] = min(self.n_multiδ[i, R_L], self.N_e - i)
        end
    end



    
    self.amplitude_accept_rates[1, R_L] /= (steps) # edge amp. update
    self.amplitude_accept_rates[2, R_L] /= (steps) # cont amp. update

    #adjust windows in edge amp. updates
    if self.amplitude_accept_rates[1, R_L] > 0.55
        self.∆n0[R_L] *= 1.25
        self.∆ε_0[R_L] *= 1.25
    elseif self.amplitude_accept_rates[1, R_L] < 0.45
        self.∆n0[R_L] /= 1.20
        self.∆ε_0[R_L] /= 1.20
    end

end


# Function to calculate bin averages and standard deviatios
function bin_averages(self::SAC, bins::Int64, std = true)

    self.χ2_res ./= bins
    self.n0_res ./= bins
    self.ωn0_res ./= bins
    self.edge_res ./= bins
    self.χ2_res[3] = sqrt(abs(self.χ2_res[3] - self.χ2_res[2]^2)) / sqrt(bins - 1)
    
    for R_L in self.R_L_index
        self.n0_res[3, R_L] = sqrt(abs(self.n0_res[3, R_L] - self.n0_res[2, R_L]^2)) / sqrt(bins - 1)
        self.ωn0_res[3, R_L] = sqrt(abs(self.ωn0_res[3, R_L] - self.ωn0_res[2, R_L]^2)) / sqrt(bins - 1)
        self.edge_res[3, R_L] = sqrt(abs(self.edge_res[3, R_L] - self.edge_res[2, R_L]^2)) / sqrt(bins - 1)
    end
end

# Function to write measured results to fule
function write_res(self::SAC, j::Int64, θ::Float64)
 
    if self.A_c > 0
        ω_cont0 = minimum(self.ω_array[self.N_e+2:end, 1])
        ω_contf = maximum(self.ω_array[self.N_e+2:end, 1])
    else
        ω_cont0 = 0.
        ω_contf = 0.
    end

    # χ2, edge, n0, etc. throughout anneal
    f = open(self.anneal_file, "a")
        println(f, j, ",", round(θ, digits=8), ",", self.χ2_min/self.N_τ, ",", self.χ2_res[2]/self.N_τ, ",", self.χ2_res[3]/self.N_τ, ",",
                   self.edge_res[2, 1], ", ", self.ωn0_res[2, 1], ",", self.n0_res[2, 1], ",", self.ε_0[1] , ",",
                   self.edge_res[3, 1], ", ", self.ωn0_res[3, 1], ",", self.n0_res[3, 1], ",",
                   -self.edge_res[2, 2], ", ", -self.ωn0_res[2, 2], ",", self.n0_res[2, 2], ",", self.ε_0[2] , ",",
                   self.edge_res[3, 2], ", ", self.ωn0_res[3, 2], ",", self.n0_res[3, 2], ",",
                   ω_cont0, ",", ω_contf)
    close(f)


    # acceptance rates
    accept_files = [self.accept_rate_file_1, self.accept_rate_file_2] #rightward and leftward edges

    for R_L in self.R_L_index

        a1 = sum(self.accept_rates[1, self.N_e+2:end, R_L]) / (self.N_ω - self.N_e) # cont single freq
        a2 = sum(self.∆ω_array[self.N_e+2:end, R_L]) / (self.N_ω - self.N_e) # cont ∆ω
        a3 = self.amplitude_accept_rates[2, R_L] #cont amplitude

        a4 = sum(self.accept_rates[1, 2:self.N_e+1, R_L]) / self.N_e # edge single freq
        # a4 = self.accept_rates[1, 2, R_L] # leading edge single freq
        a5 = sum(self.∆ω_array[2:self.N_e+1, R_L]) / self.N_e # edge ∆ω
        # a5 = self.∆ω_array[2, R_L] # leading edge ∆ω

        a6 = sum(self.accept_rates[2, 3:self.N_e - 1, R_L])/(self.N_e - 3) # edge multi freq
        a7 = sum(self.n_multiδ[3:self.N_e - 1, R_L])/(self.N_e - 3) # n_cluster

        a8 = self.amplitude_accept_rates[1, R_L] #edge amplitude
        

        f = open(accept_files[R_L], "a")
            println(f, j, ",", a1, ",", a2, ",", a3, ",", a4, ",", a5, ",", a6, ",", a7, ",", a8)
        close(f)

    end

end


# Function to write accumulated spectrum to file
function write_spec(self::SAC, n::Int64, total_steps::Int64, R_L::Int64)
    j1, j2, k = -1, -1, -1

    R_L_str = "_"*string(R_L) # _1 is leftward and _2 is rightward
  
    # standard histogram
    spec_outfile_1 = self.output_folder * "/sw" * string(n, pad=3) * R_L_str * ".dat"
    # spectral density using average delta locs and amps
    spec_outfile_2 = self.output_folder * "/dw" * string(n, pad=3) * R_L_str * ".dat"

    
    sampled_spec_hist = zeros(Float64, (self.ωi_m + 1 - self.ωi_0, 2))
    sampled_spec_bin = self.sampled_spec[:, :, R_L] ./ total_steps

    # mutlitply by normalization
    # if using ω=0 symmetric kernel, this should be multiplied by 1/2
    sampled_spec_bin .*= (self.norm * pi)
    
    # re-binning accumulated spectra in larger bin width δω_h
    δω_conversion = ceil(Int64, (self.δω_h/self.δω ))
    N_h = floor(Int64, (self.ωi_m - self.ωi_0)/δω_conversion)
    
    sampled_spec_hist[1:N_h, 1] =[sum(sampled_spec_bin[(i-1)*δω_conversion + 2:(i*δω_conversion) + 1, 1]) for i in 1:N_h]
    sampled_spec_hist[1:N_h, 2] =[sum(sampled_spec_bin[(i-1)*δω_conversion + 2:(i*δω_conversion)+ 1, 2]) for i in 1:N_h]
    sampled_spec_hist ./= self.δω_h
    
    # set cut off j1, j2 where there is no spectral weight
    for i=0:N_h
        j1=i+1
        if maximum(sampled_spec_hist[i+1, :]) > 1e-10
            j1 -= 1
            break
        end
    end
    j1 = max(1, j1)

    for i=N_h:-1:0
        j2 = i+1
        if maximum(sampled_spec_hist[i+1, :]) > 1e-10
            break
        end
    end

    
    f = open(spec_outfile_1, "w")
        println(f, "omega,S,S_edge,S_cont")
        ω = self.δω_h * (j1-1-0.5) + self.ω_0
        println(f, ω, ",", 0, ",", 0 , ",", 0)
        for i=j1:j2
            ω = self.δω_h * (i-0.5) + self.ω_0
            # println(ω)
            S1 = sampled_spec_hist[i, 1] + sampled_spec_hist[i, 2]
            println(f, ω, ",", S1, ",", sampled_spec_hist[i, 1] , ",", sampled_spec_hist[i, 2])
        end
        ω = self.δω_h * (j2+1-0.5) + self.ω_0
        println(f, ω, ",", 0, ",", 0 , ",", 0)
    close(f)


    # Averaging adjacent deltas to remove oscillations

    ω_avg_bin = self.ω_avg[:, R_L] ./ total_steps
    A_avg_bin = self.A_avg[:, R_L] ./ total_steps


    j1 = self.N_e - 1

    # calculate continuum spectrum on self-generated grid ω_avg_bin
    i1, i2 = 0, 0
    for i=1:j1
       
        # edge contribution
        sampled_spec_hist[i+1, 1] = (self.norm * pi) / (ω_avg_bin[i+2] - ω_avg_bin[i+1])
        sampled_spec_hist[i+1, 1] *= 0.5 * (A_avg_bin[i+1] + A_avg_bin[i+2])
        
        # continuum contribution
        i1 = ceil(Int64, ω_avg_bin[i+1]/self.δω) - self.ωi_0
        i2 = ceil(Int64, ω_avg_bin[i+2]/self.δω) - 1 - self.ωi_0
        sampled_spec_hist[i+1, 2] = sum(sampled_spec_bin[i1:i2, 2]) / (ω_avg_bin[i+2] - ω_avg_bin[i+1])
    end

    # use final ω_avg_bin spacing to set δω_h2 for continuum deltas that go above edge tail
    δω_h2 =  (ω_avg_bin[j1+2] - ω_avg_bin[j1+1])

    δω_conversion = floor(Int64, (δω_h2/self.δω) + 0.5)

    # Normal histogram binnings for continuum that goes above tail
    k = j1
    while true
        i1 = i2 + 1
        i2 = i1 + δω_conversion
        if i2 > (self.ωi_m - self.ωi_0)
            break
        end
        k += 1 
        sampled_spec_hist[k, 1] = 0 # no edge here
        sampled_spec_hist[k, 2] = sum(sampled_spec_bin[i1:i2, 2])/δω_h2  

        if sampled_spec_hist[k, 2] < 1.e-10
            k -= 1
            break
        end
    end

    f = open(spec_outfile_2, "w")
        println(f, "omega,S,S_edge,S_cont")
        println(f, ω_avg_bin[2], ",0,0,0")
        ω = 0.
        for i=1:j1
            ω = 0.5 * (ω_avg_bin[i+1] + ω_avg_bin[i+2])
           
            println(f, ω, ",", sampled_spec_hist[i+1, 1] + sampled_spec_hist[i+1, 2], ",",
                        sampled_spec_hist[i+1, 1], ",", sampled_spec_hist[i+1, 2])
        end
        for i=(j1+1):k
            ω += δω_h2
            println(f, ω, ",", sampled_spec_hist[i+1, 1] + sampled_spec_hist[i+1, 2], ",",
                        sampled_spec_hist[i+1, 1], ",", sampled_spec_hist[i+1, 2])
        end
        println(f,  ω + δω_h2 / 2., ",0,0,0")

    close(f)

end


###########################################################################################
## Annealing Routines
###########################################################################################

# Function to run annealing routine until χ2_target is met
# can optionally output spectrum every annealing step
function anneal(self::SAC, χ2_target::Float64, bins::Int64, spec::Int64)
    θ = self.θ_0
    i_trans = self.N_anneal * 0.2
    for i=1:self.N_anneal
        
        # ramp down updating sweeps per annealing steps
        # more steps at begining when equillibrating, then flatten out 
        # can adjust where ramp occurs (if at all) bu adjusting i_tr
        if i < i_trans 
            steps = ceil(Int64, self.anneal_steps * (1 - ((5/6) * i/i_trans)))
        else
            steps = ceil(Int64, self.anneal_steps/6)
        end
        
        run_bins(self, steps, bins, θ)
        write_res(self, i, θ)

        # write spectrum every annealing step is spec option is set
        if spec == 1
            for R_L in self.R_L_index
                write_spec(self, i, bins*steps, R_L)
            end
        end

        # save config after 5 annealing steps to 'reset' to somewhat equilibrated config
        if i == 5
            self.saved_ω_array .= self.ω_array
            self.saved_∆ω_array .= self.∆ω_array
            self.saved_n_multiδ .= self.n_multiδ
        end

        # set a target χ2 if final sampling stage
        # else check if χ2 is sufficiently close to χ2
        # else reduce θ
        if self.χ2_res[2] < χ2_target
            return θ * self.f_anneal
        elseif (self.χ2_res[2] - self.χ2_min) < (self.tol * self.N_τ) 
            return θ
        else
            θ /= self.f_anneal

            self.χ2_anneal[i] = self.χ2_res[2]
        end            

    end
    return θ
    
end 


# Basical annealing routine:
# 1) Run a full anneal (setting χ2_target == 0
#        or until χ2 is sufficiently close to χ2_min
#        or for all N_anneal steps)
# 2) Reset to equilibrated config, determine new χ2_target based on criterion with a_criterion
#        and perform second anneal until χ2 ~ χ2_target,
#        marking θ_opt where this occurs
# 3) Perform final, longer smapling using sampling_steps number sweeps at θ_opt,
#        and write spectrum/results to file.
function anneal_and_sample(self::SAC)
   
    
    # Intialize array for intial configuration for main anneal
    self.saved_ω_array = similar(self.ω_array)
    self.saved_∆ω_array = similar(self.∆ω_array)
    self.saved_n_multiδ = similar(self.n_multiδ)
   

    write_log(self, "Beginning Main Anneal.")
    anneal(self, 0., self.bins, 0)
    write_log(self, "Main Anneal Finished.")
    
    χ2_target = (self.χ2_min + self.a_criterion * sqrt(2 * self.χ2_min))
    
    n_opt = argmin(abs.(self.χ2_anneal .- χ2_target))
    θ_opt = self.θ_0 / self.f_anneal^(n_opt-1)

   
    self.ω_array .= self.saved_ω_array 
    self.∆ω_array .= self.saved_∆ω_array
    self.n_multiδ .= self.saved_n_multiδ
   
    calc_Gbar!(self, 1)
    self.χ2 = calc_χ2(self)
    self.χ2_min = self.χ2


 
    write_log(self, "Beginning Final Anneal.")   
    θ_opt = anneal(self, χ2_target, self.bins, 0)
    write_log(self, "Final Anneal Finished.")
    

    write_log(self, "Beginning Final Sampling.")
    bins = 10
    run_bins(self, self.sample_steps, bins, θ_opt)

    for R_L in self.R_L_index
        write_spec(self, 0, bins*self.sample_steps, R_L)
    end


    if self.A_c > 0
        ω_cont0 = minimum(self.ω_array[self.N_e+2:end, 1])
        ω_contf = maximum(self.ω_array[self.N_e+2:end, 1])
    else
        ω_cont0 = 0.
        ω_contf = 0.
    end
   

    f = open(self.sample_file, "a")
        println(f, 0, ",", round(θ_opt, digits=8), ",", self.χ2_min/self.N_τ, ",", self.χ2_res[2]/self.N_τ, ",",self.χ2_res[3]/self.N_τ, ",",
                   self.edge_res[2, 1], ", ", self.ωn0_res[2, 1], ",", self.n0_res[2, 1], ",", self.ε_0[1] , ",",
                   self.edge_res[3, 1], ", ", self.ωn0_res[3, 1], ",", self.n0_res[3, 1], ",",
                   -self.edge_res[2, 2], ", ", -self.ωn0_res[2, 2], ",", self.n0_res[2, 2], ",", self.ε_0[2] , ",",
                   self.edge_res[3, 2], ", ", self.ωn0_res[3, 2], ",", self.n0_res[3, 2], ",",
                   ω_cont0, ",", ω_contf)
    close(f)

    write_log(self, "Final Sampling Finished.")
end

function anneal_and_scan(self::SAC, θ_1::Float64, θ_2::Float64)

    self.N_anneal = ceil(Int64, log(self.θ_0/θ_1)/log(self.f_anneal))

   
    write_log(self, "Beginning First Anneal.")

    # Store intial configuration for main anneal
    self.saved_ω_array = similar(self.ω_array)
    self.saved_∆ω_array = similar(self.∆ω_array)
    self.saved_n_multiδ = similar(self.n_multiδ)
   

    write_log(self, "Beginning Main Anneal.")
    anneal(self, 0., self.bins ÷ 2, 0)
    write_log(self, "Main Anneal Finished.")
    
    if self.χ2_min > 2*self.N_τ
        anneal_and_scan(self, θ_1, θ_2)
        return nothing
    end

    write_log(self, "Beginning Final Sampling.")
    
    f_scan = 1.15
    N_scan = ceil(Int64, log(θ_1/θ_2)/log(f_scan))
    θ = θ_1
    for i=1:N_scan
        bins = ceil(Int64, self.bins * (1 + i/N_scan))

        run_bins(self, self.sample_steps, bins, θ)

    # for R_L in self.R_L_index
    #     write_spec(self, 0, bins*self.sample_steps, R_L)
    # end
   

        f = open(self.sample_file, "a")
            println(f, i, ",", round(θ, digits=8), ",", self.χ2_min/self.N_τ, ",", self.χ2_res[2]/self.N_τ, ",",self.χ2_res[3]/self.N_τ, ",",
                       self.edge_res[2, 1], ", ", self.ωn0_res[2, 1], ",", self.n0_res[2, 1], ",", self.ε_0[1] , ",",
                       self.edge_res[3, 1], ", ", self.ωn0_res[3, 1], ",", self.n0_res[3, 1], ",",
                       -self.edge_res[2, 2], ", ", -self.ωn0_res[2, 2], ",", self.n0_res[2, 2], ",", self.ε_0[2] , ",",
                       self.edge_res[3, 2], ", ", self.ωn0_res[3, 2], ",", self.n0_res[3, 2])
        close(f)
        θ /= f_scan

        for R_L in self.R_L_index
            write_spec(self, i, bins*self.sample_steps, R_L)
        end
    end


    write_log(self, "Final Sampling Finished.")
    return nothing
end


###########################################################################################
## Utilities 
###########################################################################################

function write_log(self, message, mode="a", output_folder=false)
    if output_folder == false
        f = open(self.output_folder * "/log.txt", mode)
            println(f, string(now(), " - ", message))
        close(f)
    else
        f = open(output_folder * "/log.txt", mode)
            println(f, string(now(), " - ", message))
        close(f)
    end
end

function write_dur(self, message)
    f = open(self.output_folder * "/log.txt", "a")
        println(f, message)
    close(f)
end


macro swap(x,y)
   quote
      local tmp = $(esc(x))
      $(esc(x)) = $(esc(y))
      $(esc(y)) = tmp
    end
end



###########################################################################################
# Function to run program
###########################################################################################
function run(A_c_in=false, A_r_in=false, p_in=false, θ_1=false, θ_2=false)

    in_file = readdlm("in_edge.in")

    
    N_e, N_c = in_file[1, :]
    ω_0, ω_m, δω_h, δω = in_file[2, :]
    p, A_c, A_r = in_file[3, :]
    θ_0, f_anneal, N_anneal, a = in_file[4, :]
    anneal_steps, sample_steps, bins = in_file[5, :]
    G_file, output_folder = in_file[6, :]
    fix_edge, kernel_type = in_file[7, 1], Symbol(in_file[7, 2])
    mode = Symbol(in_file[8, 1]) #single_edge, double_edge_in, double_edge_out, double_edge_symm

    # Check if user has submitted parameter values (for scan)
    if A_c_in != false
        A_c = A_c_in
    end
    if A_r_in != false
        A_r = A_r_in
    end

    if p_in != false
        p = p_in
    end

    # Set folder named based on mode
    if mode == :single_edge
        output_folder *= "_single/"
        A_r = 1
    elseif mode == :double_edge_in
        output_folder *= "_double_in/"
    elseif mode == :double_edge_out
        output_folder *= "_double_out/"
        ω_0 = 0
    elseif mode == :double_edge_symm
        output_folder *= "_double_symm/"
        ω_0 = 0
        A_r = 0.5
    end

    settings = []
    if fix_edge != 0
        push!(settings, "fixed") 
        ω_floor = fix_edge
        fix_edge = 1
    end


    if A_c > 0.
        if N_c == 0
            N_c = ceil(Int64, 0.5 * N_e)
            println("in_edge.in doesn't have any continuum delta's, setting to  ", N_c)
        end

        @assert A_c < 1
    end
    if A_c == 0.
        N_c = 0
       
    end


    N_ω = N_e + N_c

    push!(settings, "Nw$(string(N_e))") 

    if θ_1 != false
        push!(settings, "scan")
    end

    output_folder *= join(settings, '_')


    output_folder = output_folder * @sprintf("/Ac_%.3f", A_c)
    output_folder = output_folder * @sprintf("/p_%.3f", p)


    if mode == :double_edge_in || mode == :double_edge_out
        output_folder *= @sprintf("/Ar_%.3f", A_r)
    end
    
    if kernel_type == :bosonic
        if mode != :single_edge
            throw("Only mode that can be run with the bosonic kernel is single_edge. Exiting.")
        end
        ω_0 = 0.
    end

    c = 1 - 2*p
   
    
    mkpath(output_folder)

    cp("in_edge.in", output_folder * "/in_edge.in", force = true)
    cp(G_file, output_folder * "/t.in", force = true)

    anneal_file = output_folder*"/anneal.csv"
    accept_rate_file_1 = output_folder*"/accept_rate_1.csv"
    accept_rate_file_2 = output_folder*"/accept_rate_2.csv"
    sample_file = output_folder*"/sample.csv"
    error_file = output_folder*"/error.txt"

    println("Output folder: ", output_folder)



    
    sac = SAC(
      N_ω=N_ω, N_e = N_e, A_c=A_c, A_r=A_r, ω_0 = ω_0, ω_m=ω_m, δω=δω, δω_h=δω_h, c=[c,c],
      anneal_steps=anneal_steps, sample_steps=sample_steps, bins=bins, 
        θ_0=θ_0, N_anneal=N_anneal, f_anneal=f_anneal, a_criterion=a, G_file=G_file,
        anneal_file=anneal_file, accept_rate_file_1=accept_rate_file_1, accept_rate_file_2=accept_rate_file_2,
        sample_file=sample_file, error_file=error_file, output_folder=output_folder,
        kernel_type=kernel_type, mode=mode);

    
    sac.fix_edge = fix_edge
    if sac.fix_edge == 1
        sac.ω_floor[1] = ω_floor
        sac.ω_floor[2] = 1.5*ω_floor
    else
        sac.ω_floor[1] = ω_0
    end

    if mode == :double_edge_in
        sac.ω_floor[1] = ω_0
        sac.ω_floor[2] = -ω_m
    elseif sac.fix_edge == 1
        sac.ω_floor[1] = ω_floor
        sac.ω_floor[2] = 1.5*ω_floor
    else
        sac.ω_floor[:] .= ω_0
    end
  

    write_log(sac, "N_ω, N_e, A_r, A_c, p = " *
     string(N_ω) * ", " * string(N_e) * ", "*
      string(A_r) * ", " * string(A_c) * "," * string(p) ,
      "w")

    write_log(sac, "Beginning Initialization.")



    if mode == :single_edge || mode == :double_edge_symm
        sac.R_L_index = [1]
        write_log(sac, "Only using rightward decaying spectrum.")
    else
        sac.R_L_index = [1, 2]
        write_log(sac, "Using rightward and leftward decaying spectrum (double edge).")
    end




    write_log(sac, "Reading G_file.")
    read_G!(sac)

    write_log(sac, "Initializing Kernel.")
    init_kernal!(sac)

    write_log(sac, "Initializing struct.")
    init_struct!(sac)

    write_log(sac, "Finding Best Initial Configuration.")
    init_config_dual!(sac)
    write_log(sac, "Initialization Finished.")
    
    write_log(sac, "Running Equillibration Sweeps.")
    run_bins(sac, 5*sac.anneal_steps, sac.bins, sac.θ_0)
    write_res(sac, 0, sac.θ_0)

 
  

    if θ_1 == false
        anneal_and_sample(sac)
    else
        anneal_and_scan(sac, θ_1, θ_2)
    end

    
end


if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) > 0
        A_c, A_r, p = parse.(Float64, ARGS)
        run(A_c, A_r, p)
    else
        run()
    end
end








