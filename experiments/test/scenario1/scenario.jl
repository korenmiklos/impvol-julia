@everywhere include("../config.jl")
@everywhere import Environment, Parameters
@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium

@everywhere N, J, T, S = Environment.N, Environment.J, Environment.T, Environment.S

# parameters that govern counterfactual
@everywhere include("parameters.jl")

@time results = pmap(t -> (t, period_wrapper(A_njs, parameters, t)), 1:T)
@save "results.jld" results
