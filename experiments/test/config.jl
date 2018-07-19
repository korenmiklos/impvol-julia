#@everywhere using Logging
#@everywhere Logging.configure(level=INFO)
#@everywhere using JLD


module Environment
	include("../../calibration_params.jl")
	
	function calibrate_parameters!(parameters)
		include("../../read_data.jl")
		if !isfile("data/impvol_data.jld")
			ReadData.read_all()
		end
		data = load("../../data/impvol_data.jld")

		# mock function before implemented elsewhere

		# standard deviation for each (n,j)
		parameters[:shock_stdev] = 0.1*ones(1,N,J,1)
		# AR coefficient for each (n,j)
		parameters[:AR_decay] = 0.9*ones(1,N,J,1)

		parameters[:beta] = data[:beta]

		parameters[:gamma_jk] = CalibrateParameters.compute_gammas(parameters[:beta],data[:io_values],data[:total_output],data[:output_shares],data[:intermediate_input_shares])
		# CD case
		parameters[:nu_njt] = CalibrateParameters.compute_alphas(data[:va],parameters[:beta],parameters[:gamma_jk],parameters[:weights])
 
		parameters[:S_nt] = zeros(1,N,1,T)

		d = CalibrateParameters.expenditure_shares(data[:import_shares], parameters[:numerical_zero])

		parameters[:kappa] = CalibrateParameters.trade_costs(d, parameters[:theta], parameters[:numerical_zero])

		p_sectoral = calculate_p(data[:p_sectoral_data], data[:pwt], d, parameters[:kappa], parameter[:nu_njt], parameters[:theta])
		psi = CalibrateParameters.calculate_psi(data[:va], parameters[:weights])
		B = CalibrateParameters.calculate_B(parameters[:beta], parameters[:gamma_jk])
		xi = CalibrateParameters.calculate_xi(parameters[:theta], parameters[:eta])

		parameters[:z] = CalibrateParameters.calculate_z(p_sectoral, parameters[:beta], parameters[:gamma_jk], parameters[:kappa], psi, B, d, data[:va], xi, parameters[:theta])
		# # make Gamma more diagonal
		# gamma_jk = rand(J,J) + 1.0*eye(J)
		# parameters[:gamma_jk] = gamma_jk ./ sum(gamma_jk, 1) .* (1-beta)
	end

	export N, J, T, S, parameters

	########## environment settings
	## these are needed for data -> parameters mapping
	N = 2
	J = 5
	T = 2
	S = 100

	parameters = Dict{Symbol, Any}()
	# CES parameters
	parameters[:sigma] = 0.999
	parameters[:theta] = 4.0
	parameters[:eta] = 4.0
	parameters[:xi] = CalibrateParameters.calculate_xi(parameters[:theta], parameters[:eta])
	parameters[:N], parameters[:J], parameters[:T], parameters[:S] = N, J, T, S 

	# adaptive step size. large lambda means large steps
	parameters[:inner_step_size] = exp(-0.10*(J-1)^0.75)
	# large substitution needs more dampening
	parameters[:middle_step_size] = exp(-0.275*max(1,parameters[:sigma]))
	# any deviation from sigma=1 needs more dampening
	parameters[:outer_step_size] = exp(-0.5*abs(log(parameters[:sigma])))
	# this is log points of average input price differences
	parameters[:inner_tolerance] = 0.001
	parameters[:middle_tolerance] = 0.003
	parameters[:adjustment_tolerance] = 0.003
	parameters[:outer_tolerance] = 0.005
	parameters[:numerical_zero] = 1e-12

	########## parameters common across scenarios
	## these are function of data
	# inverse of adjustment cost, 0 if cannot readjust
	parameters[:one_over_rho] = 0.01
	parameters[:w] = [0.774074394803123; -0.201004684236153; -0.135080548288772; -0.0509519648766360]
	#parameters[:io_links] = true # and these kind of scenario parameters...
	#parameters[:china] = true

	calibrate_parameters!(parameters)
end