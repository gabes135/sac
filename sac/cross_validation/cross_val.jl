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
    output_folder *= @sprintf "/%03i%s/free" rep ab

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


    println("in_file: $G_file")
    println("out_folder: $output_folder")
   
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


function peak_cv(rep, ab)

    in_file = readdlm("in_peak.in")

    N_ω, N_p = in_file[1, :]
    A_0, ω_m, δω, δω_h = in_file[2, :]
    θ_0, f_anneal, a_criterion = in_file[3, :]
    N_anneal, anneal_steps, sample_steps = in_file[4, :]
    G_folder, output_folder = in_file[5, :]
    fix_edge, symm, kernel_type = in_file[6, 1], in_file[6, 2], Symbol(in_file[6, 3])
    
    G_file = @sprintf "%s/t_%03i%s.in" G_folder rep ab
    output_folder *= @sprintf "/%03i%s/peak" rep ab
    


    if kernel_type == :bosonic
        symm = 0

    elseif symm == 1
        output_folder *= "_symm"
    end


    if fix_edge != 0
        ω_0 = fix_edge
        fix_edge = 1
    else
        ω_0 = 0
    end

    output_folder = output_folder * @sprintf("/Np_%02i/A0_%.3f", N_p, A_0)


    mkpath(output_folder)

    cp("in_peak.in", output_folder * "/in_peak.in", force = true)
    cp(G_file, output_folder * "/t.in", force = true)

    
    println("in_file: $G_file")
    println("out_folder: $output_folder")

    sac = SAC(N_ω=N_ω, N_p = N_p, ω_m=ω_m, δω=δω, δω_h=δω_h, A_0=A_0,
            anneal_steps=anneal_steps, sample_steps=sample_steps,
            θ_0=θ_0, N_anneal=N_anneal, f_anneal=f_anneal, a_criterion=a_criterion,
            G_file=G_file,
            anneal_file=string(output_folder, "/anneal.csv"),
            sample_file=string(output_folder, "/sample.csv"),
            accept_rate_file=string(output_folder, "/accept_rate.csv"), 
            output_folder=output_folder,
            fix_edge=fix_edge, kernel_type=kernel_type, symm=symm)
    sac.ω_0 = ω_0
    sac.tol = 1e-6
    G_SAC_file = output_folder * "/GSAC.csv"
    
    open(G_SAC_file, "w") do f
    end

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
    intial_sampling(sac, θ)
    if sac.χ2_min > 1000 * sac.N_τ # Run initialization again if χ2 is too high first time 
        sac.indiv_update = true
        read_G!(sac)
        initialize!(sac)
        sac.Gbar = calc_Gbar(sac)
        sac.Gbar_new = similar(sac.Gbar)
        sac.χ2 = calc_χ2(sac)
        sac.χ2_min = sac.χ2
        sac.update_windows .= sac.ω_window / 10
        intial_sampling(sac, θ)
        # sac.indiv_update = false
    end
    write_log(sac, "Initial Sampling Finished.")

    #STEP 3
    write_log(sac, "Beginning Anneal.")


    for i=1:sac.N_anneal
        
        # steps = floor(Int64, self.anneal_steps * sqrt(i)) # Ramp up number of steps per Θ
        steps = sac.anneal_steps

        adjust_windows(sac, steps, θ)
        sample(sac, steps, θ)

        open(G_SAC_file, "a") do f
            output = vcat([round(θ, digits = 8)], sac.cov * sac.Gbar)
            writedlm(f, [output], ',')
        end



        edge_p = sac.ωi_pp * sac.δω
        A0_p = sum(sac.A_array[sac.peak_p])
        Ac_p = sum(sac.A_array[sac.cont_p])
        if sac.symm == 1 || sac.kernel_type == :bosonic
            edge_n = -edge_p
            A0_n = A0_p
            Ac_n = Ac_p
        else
            edge_n = -sac.ωi_np * sac.δω
            A0_n = sum(sac.A_array[sac.peak_n])
            Ac_n = sum(sac.A_array[sac.cont_n])
        end
        
        # Only print 4, digits
        open(sac.anneal_file, "a") do f
            writedlm(f, [[i, map(x -> round(x, digits=4), 
                          [θ, sac.χ2_min/sac.N_τ, sac.sampled_χ2/sac.N_τ,
                           edge_p, edge_n, A0_p, A0_n, Ac_p, Ac_n])...]], ',')
        end

        # Write acceptance rates
        open(sac.accept_rate_file, "a") do f
        output = cat([i], 
                     map(x -> round(x, digits=4), sac.accept_rates[[1, 2, 3, 8, 9, 10]]),
                     map(x -> round(x, digits=8), sac.update_windows[[1, 2, 8, 9]] .* sac.δω),
                     map(x -> round(x, digits=4), sac.accept_rates[[4, 5]]),
                     map(x -> round(x, digits=8), sac.update_windows[[4, 5]] .* sac.δω),
                     map(x -> round(x, digits=4), sac.accept_rates[[6, 7, 11]]),
                     map(x -> round(x, digits=8), sac.update_windows[[6, 7]].* sac.δω),
                     [round(sac.update_windows[11]* sac.δω, digits=8)],
                     dims=1)
        writedlm(f, [output], ',')
    end

    sac.χ2_anneal[i] = sac.sampled_χ2

    θ /= sac.f_anneal
    

    end

    write_log(sac, "Main Anneal Finished.")

    
    
end


function edge_cv(rep, ab)

    in_file = readdlm("in_edge.in")

    
    N_e, N_c = in_file[1, :]
    ω_0, ω_m, δω_h, δω = in_file[2, :]
    p, A_c, A_r = in_file[3, :]
    θ_0, f_anneal, N_anneal, a = in_file[4, :]
    anneal_steps, sample_steps, bins = in_file[5, :]
    G_folder, output_folder = in_file[6, :]
    fix_edge, kernel_type = in_file[7, 1], Symbol(in_file[7, 2])
    mode = Symbol(in_file[8, 1]) #single_edge, double_edge_in, double_edge_out, double_edge_symm

    G_file = @sprintf "%s/t_%03i%s.in" G_folder rep ab
    output_folder *= @sprintf "/%03i%s/edge" rep ab




    # Set folder named based on mode
    if mode == :single_edge
        A_r = 1
    elseif mode == :double_edge_out
        ω_0 = 0
    elseif mode == :double_edge_symm
        ω_0 = 0
        A_r = 0.5
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


    if fix_edge != 0
        ω_floor = fix_edge
        fix_edge = 1
        output_folder *= "_fixed"
    end

    c = 1 - 2*p
   
    
    mkpath(output_folder)

    println("in_file: $G_file")
    println("out_folder: $output_folder")

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

    sac.tol = 1e-6
    
    G_SAC_file = output_folder * "/GSAC.csv"
    
    open(G_SAC_file, "w") do f
        # header = vcat(["theta"], sac.τ)
        #writedlm(f, [header], ',')

    end

    
    sac.fix_edge = fix_edge
   

    if mode == :double_edge_in
        sac.ω_floor[1] = ω_0
        sac.ω_floor[2] = -ω_m
    elseif sac.fix_edge == 1
        sac.ω_floor .= ω_floor
    else
        sac.ω_floor .= ω_0
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


    θ = sac.θ_0
    i_trans = sac.N_anneal * 0.2
    for i=1:sac.N_anneal
        
        # ramp down updating sweeps per annealing steps
        # more steps at begining when equillibrating, then flatten out 
        # can adjust where ramp occurs (if at all) bu adjusting i_tr
        if i < i_trans 
            steps = ceil(Int64, sac.anneal_steps * (1 - ((5/6) * i/i_trans)))
        else
            steps = ceil(Int64, sac.anneal_steps/6)
        end
        
        run_bins(sac, steps, sac.bins, θ)
        write_res(sac, i, θ)


        open(G_SAC_file, "a") do f
            output = vcat([round(θ, digits = 8)], sac.cov * sac.Gbar)
            writedlm(f, [output], ',')
        end


     
        θ /= sac.f_anneal

        sac.χ2_anneal[i] = sac.χ2_res[2]
                    
    end

    write_log(sac, "Main Anneal Finished.")
    
end 




if abspath(PROGRAM_FILE) == @__FILE__ 
    param = ARGS[1]
    rep = parse(Int64, ARGS[2])

    if param == "free"
        include("../free/sac_free.jl")
        free_cv(rep, "a")
        free_cv(rep, "b")
    elseif param == "peak"
        include("../peak/sac_peak.jl")
        peak_cv(rep, "a") 
        peak_cv(rep, "b")
    elseif param == "edge"
        include("../edge/sac_edge.jl")
        edge_cv(rep, "a") 
        edge_cv(rep, "b")
    end
end
