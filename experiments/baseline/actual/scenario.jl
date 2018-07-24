# FIXME: do not parallelize FileIO
@everywhere include("../config.jl")
@everywhere using Environment

@everywhere parameters = Environment.parameters

@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium

# parameters that govern counterfactual
@everywhere include("change_parameters.jl")

@time results = pmap(t -> (t, period_wrapper(parameters, t)), 1:parameters[:T])
@save "results.jld2" results
