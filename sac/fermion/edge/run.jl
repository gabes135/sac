method = ARGS[1]

if method == "double_edge"
	include("sac_edge.jl")
	run_edge()

elseif method == "double_edge_scan"
	include("sac_double_edge.jl")
	A_c = parse.(Float64, ARGS[2])
	A_r = parse.(Float64, ARGS[3])
	run_edge(A_c, A_r, false, false)
else
	error("Invalid paramaterization.")
end
	