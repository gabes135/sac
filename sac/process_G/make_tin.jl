# Program that reads in G(τ) bins and generates t.in file for SAC

using DelimitedFiles, LinearAlgebra, Statistics, Printf

# Structure used to contain all the data for the routine
module TimeCors
Base.@kwdef mutable struct TimeCor
    
    β::Float64 # Inverse temperature
    N_b::Int64 = 0 # Number of bins, will infer from τ_grid_file and cor_file
    N_boot::Int64 # Number of bootstrap samples used to calculate covariance matrix
    
    samples::Array{Float64} = Array{Float64}(undef, 0) # Array of bootstrap samples
    
    N_τ::Int64 = 0 # Number of τ_points, will infer from τ_grid_file and cor_file
    
    τ_array::Array{Float64} = Array{Float64}(undef, 0) # Array of τ values
    G::Array{Float64} = Array{Float64}(undef, 0) # Array of mean of G(τ) values
    Gnorm::Float64 = 0.0 # Normalization G(0) + G(β) used for A(ω) in SAC
    N_τ_prime::Int64 = 0 # Number of Number of τ_points after cutting off points with too large error
    τ_ind::Array{Bool} = Array{Bool}(undef, 0) # Index of points used after cuttoff
    
    τ_grid_file::String # Input file 1: List of tau points, one time per row
    cor_file::String # Input file 2: each G(τ) bin, one G(τ) value per row
                     # Each bin should be seperated by "1"
    
    out_file::String = "t.in" # Output file read in by SAC program
     
end
end

TimeCor = TimeCors.TimeCor


# Read 
function read_data(self::TimeCors.TimeCor)
    
    
    τ_array = readdlm(self.τ_grid_file)
    N_τ = length(τ_array) - 2 # minus two to account for normalization points G(0) + G(β)

    cor_data = readdlm(self.cor_file)

    
    N_b = length(cor_data) ÷ (N_τ + 3) # plus three because of G(0),G(β), and
                                       # "1" used to seperate bins
 
    # Collect all G(τ) bins in array
    G = zeros(N_τ + 2, N_b)
    
    for k=0:(N_b - 1)
        for i=0:(N_τ+2)-1
            G[i+1,k+1] = cor_data[2 + i + (k*(N_τ+3))]
        end
    end
    
    self.τ_array = τ_array
    self.G = G 
    self.N_τ = N_τ
    self.N_b = N_b
    
end

function compute_means(self::TimeCors.TimeCor, cutoff=.2)

    G_bar = zeros(self.N_τ+2) # Average of G(τ) bins
    σ = zeros(self.N_τ+2) # Error bar on G(τ)

    samples = zeros(self.N_τ + 2,self.N_boot + 1) # Array of bootstrap samples
    
    N_τ_prime = self.N_τ
    τ_ind = trues(N_τ_prime+2) 

    bootstrap(self, 0) # Generate bootstrap samples

    for i=1:self.N_τ+2
        G_bar[i] = self.samples[i, 1]
        σ[i] = sqrt(sum((self.samples[i,2:self.N_boot+1] .- self.samples[i, 1]) .^ 2)/self.N_boot)
        if G_bar[i] < 0 || (σ[i]/G_bar[i]) > cutoff # if relative eror is above cutoff, don't include tau point
            N_τ_prime -= 1
            τ_ind[i] = false 
        end     
    end
    
    self.N_τ_prime = N_τ_prime
    self.τ_ind = τ_ind
   
    
end

# Generate bootstrap samples
function bootstrap(self, norm)
    samples = zeros(self.N_τ + 2, self.N_boot + 1)
    
    # First column is mean of all bins
    for i=1:self.N_b
        samples[:, 1] .+= self.G[:, i]
    end
    
    # The other N_boot columns are bootstrap samples (i.e. average of N_b randomly selected bins)
    for j=1:self.N_boot
        for _=1:self.N_b
            r = rand(1:self.N_b)
            samples[:, j+1] .+= self.G[:, r]
        end
    end
    
    samples ./= self.N_b
    
    # Norm of each bootstrap samples is G(0) + G(τ)
    self.Gnorm = samples[1, 1] + samples[end, 1]
    if norm == 1
        for i=0:self.N_boot
            samples[:, i+1] ./= (samples[1, i+1] + samples[end, i+1])
        end
    end
    self.samples = samples
end


# Generate covariance matrix from bootstrap samples and write t.in file
function cov_matrix(self::TimeCors.TimeCor)
    
    G_bar = zeros(self.N_τ)
    σ = zeros(self.N_τ)
    cov = zeros(self.N_τ_prime, self.N_τ_prime)
    bootstrap(self, 1)
    for i=1:self.N_τ
        G_bar[i] = self.samples[i+1, 1]
        σ[i] = sqrt(sum((self.samples[i+1,2:self.N_boot+1] .- self.samples[i+1, 1]) .^ 2)/self.N_boot)
    end
    
    G_bar = G_bar[self.τ_ind[2:end-1]]
    σ = σ[self.τ_ind[2:end-1]]
    
    compute_cov!(self, cov)
    eigen_res = eigen(cov)


    # Output file
    # First line is β, N_τ, N_boot, G(0) + G(τ)
    # Followed by four columns: τ grid, mean of G(τ) bins, error bar, and eigenvalues of covariance matrix
    open(self.out_file, "w") do f
        writedlm(f, [[self.β, self.N_τ_prime, self.N_boot, self.Gnorm]], ' ')
        writedlm(f, [self.τ_array[self.τ_ind][2:self.N_τ_prime+1] G_bar[1:end] σ[1:end] sqrt.(eigen_res.values[1:self.N_τ_prime] ./ self.N_boot)], ' ')
        for i=1:self.N_τ_prime
            println(f, i)
            eig = eigen_res.vectors[:, i]
            eig ./= norm(eig)
            writedlm(f, eig, ' ')

        end
    end 
end

function compute_cov!(self::TimeCors.TimeCor, cov::Array{Float64, 2})

    self.samples = self.samples[self.τ_ind, :]
    for j=1:self.N_τ_prime
        for i=1:self.N_τ_prime
            cov[i, j] = sum((self.samples[i+1,2:self.N_boot+1] .- self.samples[i+1,1]) .*
                            (self.samples[j+1,2:self.N_boot+1] .- self.samples[j+1,1]))
        end
    end

end


function run()

    β = 10. # Change to beta of your simulation

    N_boot = 10000 # No need to change this

    # Input files:
    cor_file = "cor.dat"
    τ_grid_file = "tgrid.dat"


    tcor = TimeCor(β=β, N_boot=N_boot, τ_grid_file=τ_grid_file, cor_file=cor_file)

    read_data(tcor)
    compute_means(tcor)
    cov_matrix(tcor)
end


run()

