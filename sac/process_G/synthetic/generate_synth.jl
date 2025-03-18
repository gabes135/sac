using DelimitedFiles, LinearAlgebra, Statistics, Printf, Dates, QuadGK

module Synths
mutable struct Synth
    spec_type::Int64 #1 Gaussians in and delta/Gaussian edge
                     #2 Edge with power-law singularity, exponential tail, and Gaussians
                     #3 Edge with power-law singularity on + and - ω-axis, exponential tail, and Gaussians
                     #5 Edge with ledge, symettric about omega = 0
                     #6 Double edge, one decaying to right and one to left with weights
                     #  A_plus and A_minus at locations ω0 and ω0_n, ω0 can be < 0
    β::Float64 # Inverse temperature
    τ_max::Float64 # Maximum time to compute G(τ)
    δτ::Float64 # Smallest τ spacing
    grid_type::Int64 # Type of τ-grid (1=linear, 2=quadratic, 3=linear+quadratic)
    M::Int64 # Target number of tau points for grid_type = 5
    τ_array::Array{Float64}
    N_τ::Int64 #Number of τ points based on τ_max, δτ, and grid_type
    
    σ::Float64 # Noise level in each generated bin as a fraction of the value G(τ=0)
    ξ::Float64 # Autocorrelation time for the noise
    N_b::Int64 # Number of bins to generate
    
    ω0::Float64 # Lower spectral edge for either delta/Gaussian or power law singularity
    ω0_n::Float64 # Negative frequency lower spectral edge for either delta/Gaussian or power law singularity
    A0::Float64 # delta/Gaussian weight at edge if spec_type = 1
                # exponent in power (0 < A0 < 1) if spec_type = 2
    ω_exp::Float64 # location of exponential fall off if spec_type = 2
                   # unused if spec_type = 1
    σ0::Float64 # Gaussian width at the edge (0 if delta) if spec_type = 1
                # decay constant in exponential tail if spec_type = 2

    A_plus::Float64 # relative weight of positive freq. singularity 
    A_minus::Float64 # relative weight of negative freq. singularity 

    
    N_G::Int64 # Number of Gaussians in continuum
    ωGs::Array{Float64} # locations of N_G continuum Gaussians 
    AGs::Array{Float64} # amplitudes of N_G continuum Gaussians
    σGs::Array{Float64} # widths of N_G continuum Gaussians
    
    G0::Array{Float64} # "clean" correlation function, without noise added
    
    spec_file::String # output file for generated S(ω)
    τ_grid_file::String # output file for list of τ points
    G_bins_file::String # output file for G(τ) bins with noise added

    kernel_type::Symbol
    
    
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
    return (exp(-ω * τ) + exp(-ω * (β - τ)))
end


function init_weights(self::Synths.Synth)
    if self.spec_type == 1
        if self.σ0 > 1e-8
            self.A0 /= (sqrt(2 * pi) * self.σ0)
        end
    end

    if self.spec_type == 1 || self.spec_type == 2 || self.spec_type == 3
        for i=1:self.N_G
            self.AGs[i] /= (sqrt(2 * pi) * self.σGs[i])
        end
    end
end

function write_spec(self::Synths.Synth)
    N_ω = 10000
    
    if self.spec_type == 1
        ωm = maximum(self.ωGs .+ (8 .* self.σGs))
        if self.A_minus > 0
            ωm_n = minimum(self.ωGs .- (8 .* self.σGs))
        end
    elseif self.spec_type == 2
        ωm = maximum(self.ωGs .+ (5 .* self.σGs))
        ωm = maximum([ωm, self.ω_exp + 1/self.σ0])
    elseif self.spec_type == 3
        ωm = maximum(self.ωGs .+ (5 .* self.σGs))
        ωm = maximum([ωm, self.ω_exp + 1/self.σ0])

        ωm_n = maximum(self.ωGs .- (5 .* self.σGs))
        ωm_n = minimum([ωm_n, -ωm])
    elseif self.spec_type == 4
        ωm = maximum(self.ωGs .+ (8 .* self.σGs))
        ωm_n = -ωm
    elseif self.spec_type == 5
        ωm = 8
        ωm_n = -8

    elseif self.spec_type == 6
        ωm = self.ω0_n
        ωm_n = self.ω0
    end
    if self.spec_type == 1 &&  self.A_minus > 0
        δω_p = (ωm - self.ω0) / N_ω
        δω_n = (self.ω0_n - ωm_n) / N_ω 

        println(ωm_n, " ", self.ω0_n, " ", self.ω0, " ", ωm)


        open(self.spec_file, "w") do f
            println(f, "omega,S")
            for i=0:N_ω-1
                ω = self.ω0_n - ((N_ω-i)*δω_n)
                s = get_spec(self, ω)
                writedlm(f, [[ω s]], ',')
            end
            println(f, self.ω0_n, ",", self.A_minus*self.A0/δω_n)
            println(f, self.ω0_n+(δω_n), ",", 0.)

            println(f, self.ω0 - (δω_p), ",", 0.)
            println(f, self.ω0, ",", self.A_plus*self.A0/δω_p)
            for i=1:N_ω
                ω = self.ω0 + (i*δω_p)
                s = get_spec(self, ω)
                writedlm(f, [[ω s]], ',')
            end
        end

    elseif self.spec_type == 1 || self.spec_type == 2

        δω = (ωm - self.ω0) / N_ω    
        println(self.ω0, " ", ωm)
        
        if self.spec_type == 1 && self.σ0 > 1e-8
            ω1 = 0.0
        else
            ω1 = self.ω0
        end
        
       
        open(self.spec_file, "w") do f
            println(f, "omega,S")
            println(f, self.ω0, ",", 0.)
            println(f, self.ω0, ",", self.A0/δω)
            for i=1:N_ω
                ω = ω1 + (i*δω)
                s = get_spec(self, ω)
                writedlm(f, [[ω s]], ',')
            end
        end
    elseif self.spec_type == 3
        δω_p = (ωm - self.ω0) / N_ω
        δω_n = (self.ω0_n - ωm_n) / N_ω 

        println(ωm_n, " ", self.ω0_n, " ", self.ω0, " ", ωm)


        open(self.spec_file, "w") do f
            println(f, "omega,S")
            for i=0:N_ω
                ω = self.ω0_n - ((N_ω-i)*δω_n)
                s = get_spec(self, ω)
                writedlm(f, [[ω s]], ',')
            end

            for i=0:N_ω
                ω = self.ω0 + (i*δω_p)
           
                s = get_spec(self, ω)
                writedlm(f, [[ω s]], ',')
            end
        end
    elseif self.spec_type > 3
        δω = (ωm - ωm_n) / N_ω
        println(ωm_n, " ", ωm, " ", δω)
        eps = 1e-3
        open(self.spec_file, "w") do f
            println(f, "omega,S")
            for i=0:N_ω-1
                ω = ωm_n + (i*δω)
                
                s = get_spec(self, ω)
                writedlm(f, [[ω s]], ',')
            end
            println(f, string(ωm) * ",0")
        end
    end

end

function get_spec(self::Synths.Synth, ω::Float64, spec_type=0)
    ε = 1e-12
    s = 0

    if spec_type == 0
        spec_type = self.spec_type
    end

    if spec_type==1 
        if (ω > self.ω0) 
            for n=1:self.N_G
                s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n]^2)))
            end
        elseif (self.σ0 > 1e-8) 
            for n=1:self.N_G
                s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n] ^2)) - (self.ω0 - ω)/self.σ0)
            end


        end
        
        # if (self.σ0 > 1e-8) 
        #    s += self.A0 * exp(-((self.ω0 - ω)^2)/(2 * (self.σ0^2)))
        # end


        if (ω < self.ω0_n) 
            for n=1:self.N_G
                s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n]^2)))
            end
        end  
        # if (self.σ0 > 1e-8) 
        #    s += self.A0 * exp(-((self.ω0_n - ω)^2)/(2 * (self.σ0^2)))
        # end


    
    elseif spec_type==2 
        if (ω > self.ω0 && ω < self.ω_exp) 
            s = (ω-self.ω0)^(-self.A0)
        elseif ω > self.ω_exp
            s = ((ω-self.ω0)^(-self.A0)) * exp(-self.σ0 * (ω - self.ω_exp)^2)

        end

        if (ω > self.ω0 - ε) 
            for n=1:self.N_G
                s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n]^2)))
            end
        end

    elseif spec_type==3 
        if (ω > self.ω0 && ω < self.ω_exp) 
            s = self.A_plus*(ω-self.ω0)^(-self.A0)
            # println([1, ω, s])
        elseif (ω < self.ω0_n && ω > -self.ω_exp) 
            s = self.A_minus*(self.ω0_n-ω)^(-self.A0)
            # println([2, ω, s, self.ω0_n-ω])
        elseif ω >= self.ω_exp
            s = self.A_plus*((ω-self.ω0)^(-self.A0)) * exp(-self.σ0 * (ω - self.ω_exp)^2)
            # println([3, ω, s])
        elseif ω <= -self.ω_exp
            s = self.A_minus*((self.ω0_n-ω)^(-self.A0)) * exp(-self.σ0 * (-ω - self.ω_exp)^2)
            # println([4, ω, s])
        end
  
        if sum(self.AGs) > 0
            if (ω > self.ω0 + ε)
                for n=1:self.N_G
                    s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n]^2)))
                end
            elseif (ω < self.ω0_n - ε) 
                for n=1:self.N_G
                    s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n]^2)))
                end
            end
        end
   

    elseif spec_type==4

        

        if (ω > self.ω0 - ε) 
            for n=1:self.N_G
                s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n]^2)))
            end
        elseif (self.σ0 > 1e-8) 
            for n=1:self.N_G
                s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n] ^2)) - (self.ω0 - ω)/self.σ0)
            end
        end
        if self.spec_type ==1
            if (self.σ0 > 1e-8) 
               s += self.A0 * exp(-((self.ω0 - ω)^2)/(2 * (self.σ0^2)))
            end
        end
    elseif spec_type == 5

        if abs(ω) < self.ω0_n
            s = 0
        elseif ω > self.ω0_n && ω< self.ω0
            s = self.A_minus
        elseif ω < -self.ω0_n && ω > -self.ω0
            s = self.A_minus

        elseif (ω < -self.ω0 && ω > -self.ω_exp) 
            s = (self.A_plus)*(-self.ω0-ω)^(-self.A0)

        elseif (ω > self.ω0 && ω < self.ω_exp) 
            s = (self.A_plus)*(ω-self.ω0)^(-self.A0)

        elseif ω >= self.ω_exp
            s = self.A_plus*((ω-self.ω0)^(-self.A0)) * exp(-self.σ0 * (ω - self.ω_exp)^2)
        elseif ω <= -self.ω_exp
            s = self.A_plus*((-self.ω0-ω)^(-self.A0)) * exp(-self.σ0 * (-ω - self.ω_exp)^2)
        end



    elseif spec_type == 6
        if ω <= self.ω0 || ω >= self.ω0_n
            s = 0
        elseif (ω > self.ω0 && ω < (self.ω0+self.ω_exp)) 
            s = self.A_plus*(ω-self.ω0)^(-self.A0)
        elseif ω >= self.ω0+self.ω_exp
            s = self.A_plus*((ω-self.ω0)^(-self.A0)) * exp(-self.σ0 * (ω - (self.ω_exp+self.ω0))^2)
        end

        if (ω < self.ω0_n && ω > (self.ω0_n-self.ω_exp)) 
            s += self.A_minus*(self.ω0_n - ω)^(-self.A0)
        elseif ω <= self.ω0_n-self.ω_exp
            s += self.A_minus*((self.ω0_n - ω)^(-self.A0)) * exp(-self.σ0 * ((self.ω0_n-self.ω_exp)-ω)^2)
        end

        #if (ω > self.ω0 - ε) 
        for n=1:self.N_G
            s += self.AGs[n] * exp(-((self.ωGs[n] - ω) ^ 2)/(2 * (self.σGs[n]^2)))
        end
        #end
    end
    return s
end

function tau_grid(self::Synths.Synth)
    n_β = floor(Int64, (self.β)/self.δτ)
    n_m = floor(Int64, (self.τ_max)/self.δτ)
    N_τ = 0
    
    τ_grid = zeros(Int64)
    
    if self.grid_type == 1
        N_τ = n_m
        τ_grid = collect(0:N_τ)
        
    elseif self.grid_type == 2
        N_τ = 0
        t1 = 0
        t2 = 0
        while true
            t2 = ((N_τ + 1)^2)  ÷ 4
            if (t2 == t1)
                t2 = t1 + 1
            end
            if t2 <= n_m
                N_τ += 1
                t1 = t2
            else
                break
            end
        end
        τ_grid = zeros(Int64,N_τ+1)
        for i=2:N_τ+1
            τ_grid[i] = ((i-1)^2) ÷ 4
            if τ_grid[i] == τ_grid[i-1]
                τ_grid[i] = τ_grid[i-1] + 1
            end
        end
        N_τ = length(τ_grid)-1
    elseif self.grid_type == 3
        N_τ = n_m
        t1 = 0
        t2 = 0
        i = 0
        while true
            
            t2 = (i+1)^2 ÷ 4
            if (t2 == t1)
                t2 = t1 + 1
            end
            if t2 > n_m && t2 <= n_β
                N_τ += 1
                t1 = t2
            elseif t2 > n_β
                break
            end
            i += 1
        end
        
        τ_grid = zeros(Int64, N_τ+1)
        for i=0:n_m
            τ_grid[i+1] = i
        end
        
        N_τ = n_m
        t1 = 0
        t2 = 0
        i = 0
        while true
           
            t2 = (i+1)^2 ÷ 4
            if (t2 == t1)
                t2 = t1 + 1
            end
            if t2 > n_m && t2 <= n_β
                N_τ += 1
                τ_grid[N_τ+1] = t2
                t1 = t2
            elseif t2 > n_β
                break
            end
            i += 1
        end
    elseif self.grid_type == 4
        n_β = floor(Int64, (self.β÷2)/self.δτ)
        n_m = floor(Int64, (self.τ_max÷2)/self.δτ)

        N_τ = n_m
        t1 = 0
        t2 = 0
        i = 0
        while true
            
            t2 = (i+1)^2 ÷ 4
            if (t2 == t1)
                t2 = t1 + 1
            end
            if t2 > n_m && t2 <= n_β
                N_τ += 1
                t1 = t2
            elseif t2 > n_β
                break
            end
            i += 1
        end
        
        τ_grid = zeros(Int64, N_τ+1)
        for i=0:n_m
            τ_grid[i+1] = i
        end
        
        N_τ = n_m
        t1 = 0
        t2 = 0
        i = 0
        while true
           
            t2 = (i+1)^2 ÷ 4
            if (t2 == t1)
                t2 = t1 + 1
            end
            if t2 > n_m && t2 <= n_β
                N_τ += 1
                τ_grid[N_τ+1] = t2
                t1 = t2
            elseif t2 > n_β
                break
            end
            i += 1
        end

        τ_grid = vcat(τ_grid, reverse(1+2*self.τ_max ÷ self.δτ .- τ_grid ))
        N_τ = length(τ_grid) - 2
    elseif self.grid_type == 5
        #log spaced grid to beta//2, assuming G(beta) = G(0) so norm is 2*G(0)
        N = (self.τ_max ÷ self.δτ) + 1
        M = self.M#40
        τ_grid = Integer.(unique(round.(10 .^ range( log10(1), log10(N), length = M ), digits=0))) .- 1
        N_τ = length(τ_grid)-1
    elseif self.grid_type == 6
        N = (self.τ_max ÷ self.δτ) + 1
        M = self.M#40
        τ_grid = Integer.(unique(round.(10 .^ range( log10(1), log10(N), length = M ), digits=0))) .- 1
        τ_grid = vcat(τ_grid, reverse(1+2*self.τ_max ÷ self.δτ .- τ_grid ))
        N_τ = length(τ_grid)-1

     elseif self.grid_type == 7
        N = (self.τ_max ÷ self.δτ) + 1
        M = self.M#40
        τ_grid = Integer.(unique(round.(10 .^ range( log10(1), log10(N), length = M ), digits=0))) .- 1
        τ_grid = vcat(τ_grid, (2*n_β) .- reverse(τ_grid) .- 1, [self.β])
        N_τ = length(τ_grid)-1
    end
    
    println("τ_max = ", round(τ_grid[end]* self.δτ, digits =2))
    self.τ_array = τ_grid * self.δτ
    self.N_τ = N_τ

    
end

function integrand(self::Synths.Synth, τ, ω, spec_type = 0)

    if self.kernel_type == :finiteT
        K_function = finiteT_K
    elseif self.kernel_type == :zeroT
        K_function = zeroT_K
    elseif self.kernel_type == :bosonic
        K_function = bosonic_K
    else
        throw("Invalid Kernel type.")
    end

    return get_spec(self, ω, spec_type) * K_function(ω, τ, self.β)

end



function make_G_tau(self::Synths.Synth)

    if self.kernel_type == :finiteT
        K_function = finiteT_K
    elseif self.kernel_type == :zeroT
        K_function = zeroT_K
    elseif self.kernel_type == :bosonic
        K_function = bosonic_K
    else
        throw("Invalid Kernel type.")
    end
    

    ωm = 0.0
    ωm_n = 0.0
    if self.spec_type == 1 
        ωm = maximum(self.ωGs .+ (10 .* self.σGs))
        if self.A_minus > 0
            ωm_n = minimum(self.ωGs .- (10 .* self.σGs))
        end
    elseif self.spec_type == 2 
        ωm = maximum(self.ωGs .+ (10 .* self.σGs))
        ωm = maximum([ωm, self.ω_exp + 10/self.σ0])
    elseif self.spec_type == 3 
        ωm = maximum(self.ωGs .+ (5 .* self.σGs))
        ωm = maximum([ωm, self.ω_exp + 0.75/self.σ0])
        ωm_n = -ωm

    elseif self.spec_type == 4
        ωm = maximum(self.ωGs .+ (10 .* self.σGs))
        ωm_n = -ωm
    elseif self.spec_type == 5
        ωm = 8
        ωm_n = -ωm
    elseif self.spec_type == 6
        ωm_n = self.ω0
        ωm = self.ω0_n
    end

    
    G0 = zeros(self.N_τ+1)
    
    
    for j=0:self.N_τ
        τ = self.τ_array[j+1]
        
        #println(τ)
        G0_val = 0.0
        if self.spec_type==1 
            if self.σ0 < 1e-8
                G0_val = quadgk(ω -> integrand(self, τ, ω), self.ω0, ωm, rtol=1e-13)[1]
                
                if self.A_minus > 0 
                    G0_val += quadgk(ω -> integrand(self, τ, ω), ωm_n, self.ω0_n, rtol=1e-13)[1]
                end
            else
                G0_val = quadgk(ω -> integrand(self, τ, ω), ωm_n, ωm, rtol=1e-13)[1]
            end
            
            if self.σ0 < 1e-8
                G0_val += self.A_plus*self.A0* K_function(self.ω0, τ, self.β)
                G0_val += self.A_minus*self.A0* K_function(self.ω0_n, τ, self.β)
            end
        elseif self.spec_type==2
            ω1 = 0.05
            # integrate power-law part for ω ~ ω0 close to edge analytically
            G0_val += small_ω(self, τ, ω1) * exp(-τ*self.ω0) # dont forget the kernal

            ω1 += self.ω0
            # integrate rest of spectrum
            G0_val += quadgk(ω -> integrand(self, τ, ω), ω1, ωm, rtol=1e-13)[1] 
            # integrate Guassian part for ω ~ ω0 close to edge
            G0_val += quadgk(ω -> integrand(self, τ, ω, 1), self.ω0, ω1, rtol=1e-13)[1] 
            
        
        elseif  self.spec_type == 3 || self.spec_type == 4
            a = quadgk(ω -> integrand(self, τ, ω), ωm_n, -self.ω0, rtol=1e-13)[1]
            b = quadgk(ω -> integrand(self, τ, ω), self.ω0, ωm, rtol=1e-13)[1]
            # println([a, b])
            G0_val = a#quadgk(ω -> integrand(self, τ, ω), ωm_n, -self.ω0, rtol=1e-13)[1]
            G0_val += b#quadgk(ω -> integrand(self, τ, ω), self.ω0, ωm, rtol=1e-13)[1]

        elseif  self.spec_type == 5
            a = quadgk(ω -> integrand(self, τ, ω), ωm_n, -self.ω0, rtol=1e-13)[1]
            b = quadgk(ω -> integrand(self, τ, ω), -self.ω0, 0, rtol=1e-13)[1]
            c = quadgk(ω -> integrand(self, τ, ω), 0, self.ω0, rtol=1e-13)[1]
            d = quadgk(ω -> integrand(self, τ, ω), self.ω0, ωm, rtol=1e-13)[1]
            # println([a, b])
            G0_val = a+b+c+d#quadgk(ω -> integrand(self, τ, ω), ωm_n, -self.ω0, rtol=1e-13)[1]
        elseif  self.spec_type == 6
            G0_val = quadgk(ω -> integrand(self, τ, ω), ωm_n, ωm, rtol=1e-13)[1]
        end


        G0[j+1] = G0_val/pi

    end
    
    self.G0 = G0


    # remove_inds = []
    # for i=1:length(self.τ_array)
    #     if G0[i] < synth.σ
    #         push!(remove_inds, i)
    #     end
    # end
    # self.N_τ -= length(remove_inds)
    # deleteat!(self.τ_array, remove_inds)
    # deleteat!(self.G0, remove_inds)

    open(self.τ_grid_file, "w") do f
        # writedlm(f, self.N_τ, ' ')
        writedlm(f, round.(synth.τ_array, digits=8), ' ')
    end
    #for i=1:length(remove_inds)
    

end

function add_noise(self::Synths.Synth)
    if self.spec_type == 1 && self.A_minus > 0
        σ = self.σ * (self.G0[1] + self.G0[end])
    elseif self.spec_type == 4
        σ = self.σ * (self.G0[1] + self.G0[end])
    elseif self.spec_type == 6 && self.ω0 < 0
        σ = self.σ * self.G0
    elseif self.spec_type == 6 && self.ω0 > 0
        σ = self.σ * self.G0[1]
    else
        σ = self.σ * (self.G0[1] * 2)
    end

    σG = σ .* randn(self.N_τ+1) 
    #σG = σ * randn() 
    
    σG_corr = zeros(self.N_τ+1)
    
    for j=0:self.N_τ
        σ_corr_norm = 0
        for i=0:self.N_τ
            exp_factor = exp(-abs(self.τ_array[i+1] - self.τ_array[j+1])/self.ξ)
            σG_corr[j+1] += (exp_factor * σG[i+1])
            σ_corr_norm += exp_factor^2
        end
        σG_corr[j+1] /= sqrt(σ_corr_norm)
    end
    return self.G0 .+ σG_corr
end

function write_Gbins(self)
    open(self.G_bins_file, "w") do f
        for j=1:self.N_b
            G1 = add_noise(self)
            writedlm(f, "1", ' ')
            writedlm(f, round.(G1, digits=8), ' ')
        end
    end
end
    


function Synth(spec_type, β, τ_max, δτ, grid_type, M,
               σ, ξ, N_b,
               ω0, ω0_n, A0, ω_exp, σ0,
               A_plus, A_minus,
               N_G, ωGs, AGs, σGs,
               spec_file, τ_grid_file, G_bins_file, kernel_type)

    τ_array = zeros()
    N_τ = 0
    G0 = zeros()


    args = [spec_type, β, τ_max, δτ, grid_type, M, τ_array, N_τ,
     σ, ξ, N_b,
    ω0, ω0_n, A0, ω_exp, σ0,
    A_plus, A_minus,
    N_G, ωGs, AGs, σGs,
    G0, spec_file, τ_grid_file, G_bins_file, kernel_type]

    return Synths.Synth(args...)

end

##################################################################

function small_ω(self, τ, ω)
    res = (ω^(1-self.A0))/(1-self.A0)
    f1 = 1
    for n=1:200
        f1*=(-τ/n)
        f2 = (f1/((n+1) - self.A0))*(ω^((n+1)-self.A0))
        res += f2
        if abs(f2) < 1e-20
            break
        end
    end
    return res
end

