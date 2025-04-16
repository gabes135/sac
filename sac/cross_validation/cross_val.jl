using Dates, DelimitedFiles, LinearAlgebra, Printf




function free_cv(rep, ab)

    in_file = readdlm("in_free.in")

    par = in_file[1,1]
    N_ω, ω_0, ω_m, δω, δω_h = in_file[2, :]
    θ_0, f_anneal, f_final, a1, a2 = in_file[3, :]
    N_anneal, anneal_steps, sample_steps = in_file[4, :]
    G_folder, output_folder = in_file[5, :]
    symm, kernel_type = in_file[6, 1], Symbol(in_file[6, 2])

    G_file = @sprintf "%s/t_%03i%s.in" G_folder rep ab
    output_folder *= @sprintf "/%03i%s" rep ab

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

    sac.tol = 1e-6
    
    G_SAC_file = output_folder * "/GSAC.csv"
    
    open(G_SAC_file, "w") do f
        # header = vcat(["theta"], sac.τ)
        #writedlm(f, [header], ',')

    end
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



    θ = sac.θ_0
    for i=1:sac.N_anneal
        steps = sac.anneal_steps
        
        adjust_windows(sac, steps, θ)
        sample(sac, steps, θ)

        open(G_SAC_file, "a") do f
            output = vcat([round(θ, digits = 8)], sac.cov * sac.Gbar)
            writedlm(f, [output], ',')
        end

    
        open(sac.anneal_file, "a") do f
           writedlm(f, [[i, round(θ, digits = 8),
                            round(sac.χ2_min/sac.N_τ, digits = 4),
                            round(sac.sampled_χ2/sac.N_τ, digits = 4)]], ',')
        end
        
        open(sac.accept_rate_file, "a") do f
            rounded_values = round.(vcat(sac.accept_rates..., sac.update_windows[[1, 2, 4]] .* sac.δω), digits=8)

            output = cat([i], rounded_values, dims=1)

            writedlm(f, [output], ',')
        end

        sac.χ2_anneal[i] = sac.sampled_χ2
        
        # if wr is set to true then spectrum is output at each temperature in anneal
        # instead of just during final sampling step 
        
        # write_spec(sac, i)
       

        # Exit main anneal is sufficiently close to χ2_min
        if (sac.sampled_χ2 - sac.χ2_min) < (sac.tol * sac.N_τ) 
            return
        else
            θ /= sac.f_anneal
        end  

    end

    write_log(sac, "Main Anneal Finished.")
    
end 




if abspath(PROGRAM_FILE) == @__FILE__ 
    param = ARGS[1]
    rep = parse(Int64, ARGS[2])

    if param == "free"
        include("sac_free.jl")
        free_cv(rep, "a")
        free_cv(rep, "b")
    elseif param == "peak"
        include("sac_peak.jl")
    elseif param == "edge"
        include("sac_edge.jl") 
    end
end
