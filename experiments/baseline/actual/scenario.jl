@everywhere include("../config.jl")
@everywhere using Environment
@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium

# parameters that govern counterfactual
@everywhere include("change_parameters.jl")

@time results = pmap(t -> (t, period_wrapper(draw_next_productivity(parameters,t), parameters, t)), 1:parameters[:T])
@save "results.jld" results
