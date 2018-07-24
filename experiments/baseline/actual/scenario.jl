# FIXME: do not parallelize FileIO
include("../config.jl")
using Environment

include("../../../equilibrium.jl")
using ImpvolEquilibrium

# parameters that govern counterfactual
include("change_parameters.jl")

@time results = pmap(t -> (t, period_wrapper(draw_next_productivity(parameters,t), parameters, t)), 1:parameters[:T])
@save "results.jld2" results
