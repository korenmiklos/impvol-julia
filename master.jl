# set random seed for reproducability
srand(2311)

@everywhere include("utils.jl")
@everywhere include("equilibrium.jl")
@everywhere using ImpvolEquilibrium

# for matrix conformity, store all variables in a 4-dimensional array:
# mnjt: destination, source, sector, time

# per-period random variables are stored as
# mnjs: destination, source, sector, state

# load parameters from config file
@everywhere include("experiments/test/config.jl")
@everywhere include("experiments/test/scenario1/parameters.jl")
@everywhere import Environment, Parameters, Scenario

parameters = merge(Environment.parameters, Parameters.parameters, Scenario.parameters)
coerce_parameters!(parameters)

A_njs = 1.0 .+ rand(1,Environment.N,Environment.J,Environment.S)
# to test CES
A_njs[1,1,1,1] = 1.0
A_njs[1,1,2,1] = 0.5

@time results = pmap(t -> (t, period_wrapper(A_njs, parameters, t)), 1:Environment.T)

using JLD
@save "experiments/test/scenario1/results.jld" results
