using Dates, DelimitedFiles, LinearAlgebra, Printf

function read_G!(G_file)
   
    t_in = readdlm(G_file)
    
    β, N_τ, n_bins, norm = t_in[1,:]
    
    N_τ = Int(N_τ)
    
    cov = zeros((N_τ, N_τ))
    
    τ, G, σ = t_in[2:N_τ+1,1], t_in[2:N_τ+1,2], t_in[2:N_τ+1,4]
    σ_inv = 1. ./ σ
    
    cov_start = N_τ+2
    
    for j=0:N_τ-1
        @assert j+1 == t_in[cov_start + j*(N_τ+1),1]
        cov[:, j+1] = t_in[cov_start + j*(N_τ+1) + 1:cov_start + j*(N_τ+1) + N_τ,1]
    end
    
 
    G_D = transpose(cov) * G 
    
    return G_D, cov, σ_inv, N_τ
end 

function read_GSAC(GSAC_file, G_D, cov, σ_inv, N_τ)
    GSAC = readdlm(GSAC_file, ',')
    θs, GSAC = GSAC[:, 1], GSAC[:, 2:end]

    println(size(GSAC))
    N_θ = length(θs)
  
    χ2s = []
    for t=1:N_θ
        GSAC_t = transpose(cov) * GSAC[t, :]
        push!(χ2s, calc_χ2(GSAC_t, G_D, σ_inv, N_τ))
    end

    return χ2s

end

function calc_χ2(GSAC, G_D, σ_inv, N_τ)
    χ2 = 0
    @inbounds @simd for i=1:N_τ
        χ2 += ((G_D[i] - GSAC[i]) * σ_inv[i]) ^ 2
    end
    return χ2
end


function run(G_file, GSAC_file)
    G_D, cov, σ_inv, N_τ = read_G!(G_file)

    χ2s = read_GSAC(GSAC_file, G_D, cov, σ_inv, N_τ) ./ N_τ
end


if abspath(PROGRAM_FILE) == @__FILE__
    G_file = "in_files/hchain_beta2048/t_001a.in"
    GSAC_file = "out_files/hchain_beta2048/001a/GSAC.csv"
    χ2s = run(G_file, GSAC_file)
    println(χ2s)
end





