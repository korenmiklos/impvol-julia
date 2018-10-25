# load ImpvolEquilibrium first so that methods are accessible
@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium

@everywhere include("../init_parameters.jl")

# parameters that govern counterfactual
@everywhere include("change_parameters.jl")


@time results = pmap(t -> (t, ImpvolEquilibrium.period_wrapper(parameters, t)), 1:parameters[:T])
CalibrateParameters.jld_saver(results)
