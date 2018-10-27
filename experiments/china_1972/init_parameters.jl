	## these are needed for data -> parameters mapping
	parameters = Dict{Symbol, Any}()

	# CES parameters
	parameters[:sigma] = 0.999
	parameters[:theta] = 4.0
	parameters[:eta] = 4.0

	########## parameters common across scenarios
	## these are function of data
	# inverse of adjustment cost, 0 if cannot readjust
	parameters[:one_over_rho] = 0.001

	include("../config.jl")
	# change parameters after reading data, but common across scenarios
	for t = 2:parameters[:T]
		parameters[:kappa_mnjt][5,:,:,t] = parameters[:kappa_mnjt][5,:,:,1]
		parameters[:kappa_mnjt][:,5,:,t] = parameters[:kappa_mnjt][:,5,:,1]
	end
