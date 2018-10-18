	## these are needed for data -> parameters mapping
	parameters = Dict{Symbol, Any}()

	# CES parameters
	parameters[:sigma] = 0.999
	parameters[:theta] = 4.0
	parameters[:eta] = 4.0

	########## parameters common across scenarios
	## these are function of data
	# inverse of adjustment cost, 0 if cannot readjust
	parameters[:one_over_rho] = 0.01

	include("../config.jl")
	# change parameters after reading data, but common across scenarios
	## Balanced trade
	parameters[:beta_j] = ones(size(parameters[:beta_j]))
	parameters[:gamma_jk] = zeros(size(parameters[:gamma_jk]))
	using .CalibrateParameters
	parameters[:B] = CalibrateParameters.calculate_B(parameters)
	using FileIO
	data = load("../../../data/impvol_data.jld2")
	parameters[:nu_njt] = CalibrateParameters.compute_alpha(parameters, data)
