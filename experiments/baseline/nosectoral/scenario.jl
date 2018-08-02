@everywhere include("../init_parameters.jl")
@everywhere using Environment
@everywhere parameters = Environment.parameters

# parameters that govern counterfactual
@everywhere include("change_parameters.jl")

@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium

@time results = pmap(t -> (t, period_wrapper(parameters, t)), 1:parameters[:T])
Environment.CalibrateParameters.jld_saver(results)