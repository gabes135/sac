include("base.jl")
using .SACs


mode = ARGS[1]

if mode == "free"
    include("sac_free.jl")
    using .FreeSACs
    FreeSACs.run()
elseif mode == "peak"
    include("sac_peak.jl")
    using .PeakSACs
    PeakSACs.run()
elseif mode =="edge"
    include("sac_edge.jl")
    using .EdgeSACs
    EdgeSACs.run()
else
    error("Invalid mode! Please select 'free', 'peak', or 'edge'.")

end


# if abspath(PROGRAM_FILE) == @__FILE__
    
# end