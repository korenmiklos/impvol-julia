module Environment
	export N, J, T, S, parameters
	########## environment parameters
	N = 2
	J = 5
	T = 2
	S = 100

	parameters = Dict{Symbol, Any}()
	# CES parameters
	## FIXME: nu has to be calibrated to data
	parameters[:nu_njt] = ones(1, N, J, T)
	parameters[:sigma] = 0.999
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
	parameters[:numerical_zero] = 1e-6
end

@everywhere using Logging
@everywhere Logging.configure(level=INFO)

module Parameters
	import ..Environment: N, J, T
	export parameters
	parameters = Dict{Symbol, Any}()
	########## parameters common across scenarios
	# inverse of adjustment cost, 0 if cannot readjust
	parameters[:one_over_rho] = 0.01

	# standard deviation for each (n,j)
	parameters[:shock_stdev] = 0.1*ones(1,N,J,1)
	# AR coefficient for each (n,j)
	parameters[:AR_decay] = 0.9*ones(1,N,J,1)

	beta = 0.25 + 0.75 * rand(1, J)

	
	parameters[:beta]=beta
	parameters[:theta]=4
	parameters[:S_nt]=zeros(1,N,1,T)

	# make Gamma more diagonal
	gamma_jk = rand(J,J) + 1.0*eye(J)
	parameters[:gamma_jk] = gamma_jk ./ sum(gamma_jk, 1) .* (1-beta)
end