# load ImpvolEquilibrium first so that methods are accessible
@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium

@everywhere include("../init_parameters.jl")

# parameters that govern counterfactual
@everywhere include("change_parameters.jl")
@everywhere parameters[:S] = 2

@everywhere srand(7094)

@time results = pmap(t -> (t, ImpvolEquilibrium.period_wrapper(parameters, t)), 1:2)
CalibrateParameters.jld_saver(results)
