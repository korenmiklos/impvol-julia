# load ImpvolEquilibrium first so that methods are accessible
include("../../equilibrium.jl")
using Logging, ImpvolEquilibrium, FileIO, JLD2

## these are needed for data -> parameters mapping
parameters = Dict{Symbol, Any}()

# CES parameters
parameters[:sigma] = 0.999
parameters[:theta] = 4.0
parameters[:eta] = 4.0

########## parameters common across scenarios
## these are function of data
# inverse of adjustment cost, 0 if cannot readjust
parameters[:one_over_rho] = 0.0

Logging.configure(level=DEBUG)

# save counterfactual trade shares
data = load("../../data/impvol_data.jld2")
data["import_shares"] = data["import_shares"][:,:,:,1:1] .* ones(data["import_shares"])
save("../../data/data_fixed_trade.jld2", data)

include("../../calibrate_params.jl")
using CalibrateParameters

if !haskey(parameters, :S)
	parameters[:S] = 101
end

parameters[:numerical_zero] = 1e-12

parameters[:bp_weights] = [0.774074394803123; -0.201004684236153; -0.135080548288772; -0.0509519648766360]

# set random seed so that all scenarios comparable
srand(9181)
CalibrateParameters.calibrate_parameters!(parameters, "../../data/data_fixed_trade.jld2")

# adaptive step size. large lambda means large steps
parameters[:inner_step_size] = exp(-0.06*(parameters[:J]-1)^0.75)
# large substitution needs more dampening
parameters[:middle_step_size] = exp(-0.8*max(1,parameters[:sigma]))
parameters[:adjustment_step_size] = 0.25
# any deviation from sigma=1 needs more dampening
parameters[:outer_step_size] = exp(-0.5*(1+abs(log(parameters[:sigma]))))
# this is log points of average input price differences
parameters[:inner_tolerance]  = 0.0002
parameters[:middle_tolerance] = 0.0003
parameters[:adjustment_tolerance] = 0.0004
parameters[:outer_tolerance] = 0.001

# maximum number of iterations in each loop
parameters[:max_iter_inner] = 200
parameters[:max_iter_middle] = 100
parameters[:max_iter_adjustment] = 100
parameters[:max_iter_outer] = 100

# change parameters after reading data, but common across scenarios
jldopen("common_parameters.jld2", "w") do file
	file["parameters"] = parameters
end
