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
    
    folder = ARGS[1]
    reps = parse(Int64, ARGS[2])
    param = ARGS[3]

    out_folder = @sprintf "out_files/%s/chi2" folder 
    mkpath(out_folder)

    anneal_file = @sprintf "out_files/%s/%03ib/%s/anneal.csv" folder 1 param
    theta = readdlm(anneal_file, ',', skipstart=1)[:, 2]

    out_file_s = @sprintf "%s/%s_s.csv" out_folder split(param, "/")[1]
    open(out_file_s, "w") do f
        writedlm(f, [theta], ',')
    end

    out_file_v = @sprintf "%s/%s_v.csv" out_folder split(param, "/")[1]
    open(out_file_v, "w") do f
        writedlm(f, [theta], ',')
    end
    
    for rep in 1:reps
        G_file = @sprintf "in_files/%s/t_%03ia.in" folder rep
        GSAC_file = @sprintf "out_files/%s/%03ib/%s/GSAC.csv" folder rep param
        χ2_v_ab = run(G_file, GSAC_file)

        χ2_s_file = @sprintf "out_files/%s/%03ib/%s/anneal.csv" folder rep param
        χ2_s_b = readdlm(χ2_s_file, ',', skipstart=1)[:, 4]


        G_file = @sprintf "in_files/%s/t_%03ib.in" folder rep
        GSAC_file = @sprintf "out_files/%s/%03ia/%s/GSAC.csv" folder rep param
        χ2_v_ba = run(G_file, GSAC_file)

        χ2_s_file = @sprintf "out_files/%s/%03ia/%s/anneal.csv" folder rep param
        χ2_s_a = readdlm(χ2_s_file, ',', skipstart=1)[:, 4]

        χ2_v = (χ2_v_ab .+ χ2_v_ba) ./ 2
        χ2_s = (χ2_s_a .+ χ2_s_b) ./ 2

        open(out_file_v, "a") do f
            writedlm(f, [χ2_v], ',')
        end
        open(out_file_s, "a") do f
            writedlm(f, [χ2_s], ',')
        end
    end
    
end




