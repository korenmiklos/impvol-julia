module Environment
	using Logging
	Logging.configure(level=DEBUG)
	include("../../calibrate_params.jl")
	using .CalibrateParameters
	
	export parameters

	########## environment settings
	## these are needed for data -> parameters mapping
	parameters = Dict{Symbol, Any}()
	parameters[:S] = 100

	# CES parameters
	parameters[:sigma] = 0.999
	parameters[:theta] = 4.0
	parameters[:eta] = 4.0

	parameters[:numerical_zero] = 1e-12

	########## parameters common across scenarios
	## these are function of data
	# inverse of adjustment cost, 0 if cannot readjust
	parameters[:one_over_rho] = 0.01
	parameters[:bp_weights] = [0.774074394803123; -0.201004684236153; -0.135080548288772; -0.0509519648766360]
	#parameters[:io_links] = true # and these kind of scenario parameters...
	#parameters[:china] = true

	CalibrateParameters.calibrate_parameters!(parameters)

	# adaptive step size. large lambda means large steps
	parameters[:inner_step_size] = exp(-0.10*(parameters[:J]-1)^0.75)
	# large substitution needs more dampening
	parameters[:middle_step_size] = exp(-0.275*max(1,parameters[:sigma]))
	# any deviation from sigma=1 needs more dampening
	parameters[:outer_step_size] = exp(-0.5*abs(log(parameters[:sigma])))
	# this is log points of average input price differences
	parameters[:inner_tolerance] = 0.001
	parameters[:middle_tolerance] = 0.003
	parameters[:adjustment_tolerance] = 0.003
	parameters[:outer_tolerance] = 0.005
end