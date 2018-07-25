
module CalibrateParameters
	include("calibration_utils.jl")
	include("utils.jl")
	using FileIO
	
	function calibrate_parameters!(parameters)
		data = load("../../../data/impvol_data.jld2")

		_, N, J, T = size(data["beta"])
		parameters[:N], parameters[:J], parameters[:T] = N, J, T
		parameters[:beta_j] = mean(data["beta"],(1,2,4))

		# standard deviation for each (n,j)
		parameters[:shock_stdev] = 0.1*ones(1,N,J,1)
		# AR coefficient for each (n,j)
		parameters[:AR_decay] = 0.9*ones(1,N,J,1)

		parameters[:gamma_jk] = compute_gamma(parameters, data)
		# CD case
		parameters[:nu_njt] = compute_alpha(parameters, data)
 
		parameters[:S_nt] = zeros(1,N,1,T)

		parameters[:d] = expenditure_shares(parameters, data)

		parameters[:kappa_mnjt] = trade_costs(parameters)

		parameters[:p_sectoral] = calculate_p(parameters, data)
		parameters[:psi] = calculate_psi(parameters, data)
		parameters[:B_j] = calculate_B(parameters)
		parameters[:xi] = calculate_xi(parameters)

		parameters[:z] = calculate_z(parameters, data)

		parameters[:A] = calculate_A(parameters)

		parameters[:A_njs] = map(t -> (t, draw_next_productivity(parameters, t)), 1:parameters[:T])
	end

	function compute_gamma(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		beta = parameters[:beta_j]
		io_values = data["io_values"]
		total_output = data["total_output"]
		output_shares = data["output_shares"]
		intermediate_input_shares = data["intermediate_input_shares"]

		# Summing sectors
		# Agriculture & mining
		io_values[1,:,:,:]  = sum(io_values[1:2,:,:,:],1)
		io_values         = io_values[setdiff(1:end,2),:,:,:]
		io_values[:,1,:,:]  = sum(io_values[:,1:2,:,:],2)
		io_values         = io_values[:,setdiff(1:end,2),:,:]
		total_output[:,1,:,:] = sum(total_output[:,1:2,:,:],2)
		total_output      = total_output[:,setdiff(1:end,2),:,:]

		# Services
		io_values[18,:,:,:]  = sum(io_values[18:end,:,:,:],1)
		io_values          = io_values[setdiff(1:end,19:end),:,:,:]
		io_values[:,18,:,:]  = sum(io_values[:,18:end,:,:],2)
		io_values          = io_values[:,setdiff(1:end,19:end),:,:]
		total_output[:,18,:,:] = sum(total_output[:,18:end,:,:],2)
		total_output       = total_output[:,setdiff(1:end,19:end),:,:]

		# Split rows
		dupl_idx = [1, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 14, 15, 16, 17, 18]
		io_values_dupl = io_values[dupl_idx,:,:,:]

		output_shares_full = ones(size(io_values_dupl))
		split_idx = [2, 3, 4, 5, 6, 8, 9, 17, 18, 19]
		output_shares_full[split_idx,:,:,:] = permutedims(repeat(output_shares, outer = [1,1,size(io_values_dupl,2),1]),(2,3,1,4))

		io_values_new = io_values_dupl .* output_shares_full

		# Split columns
		io_values_dupl = io_values_new[:,dupl_idx,:,:]

		intermediate_input_shares_full = ones(size(io_values_dupl))
		intermediate_input_shares_full[:,split_idx,:,:] = repeat(intermediate_input_shares, outer = [size(io_values_dupl,1),1,1,1])

		io_values_new = io_values_dupl .* intermediate_input_shares_full

		total_output = total_output[:,dupl_idx,:,:]
		output_shares_full = ones(size(total_output))
		output_shares_full[:,split_idx,:,:] = output_shares
		total_output = total_output .* output_shares_full

		# Correct order of sectors
		io_values_new[18,:,:,:], io_values_new[19,:,:,:], io_values_new[20,:,:,:] = io_values_new[20,:,:,:], io_values_new[18,:,:,:], io_values_new[19,:,:,:]
		io_values_new[:,18,:,:], io_values_new[:,19,:,:], io_values_new[:,20,:,:] = io_values_new[:,20,:,:], io_values_new[:,18,:,:], io_values_new[:,19,:,:]
		total_output[:,18,:,:], total_output[:,19,:,:], total_output[:,20,:,:] = total_output[:,20,:,:], total_output[:,18,:,:], total_output[:,19,:,:]

		# Compute gamma
		gamma = io_values_new ./ repeat(total_output, outer = [size(io_values_new,1),1,1,1])
		gamma = mean(gamma,4)
		gamma = gamma .* permutedims(1-beta,(1,3,2,4)) ./ sum(gamma,1)
		return gamma = squeeze(gamma,(3,4))
	end

	function compute_alpha(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		va = data["va"]
		beta = parameters[:beta_j]
		gamma = parameters[:gamma_jk]
		weights = parameters[:bp_weights]

		alpha = zeros(J,T)

		for t in 1:T
			va_t = transpose(sum(va[1,:,:,t],1))
			alpha[:,t] = (eye(J) - gamma) * diagm(1 ./ beta[:],0) * va_t / sum(va_t)
		end

		# Replace negative elements with 0
		alpha = (alpha + abs.(alpha)) / 2

		# Smooth the series
		alpha_c, alpha_t = DetrendUtilities.detrend(alpha, weights)

		# Normalization
		alpha = alpha_t ./ sum(alpha_t,1)

		return alpha = permutedims(cat(ndims(alpha) + 2,alpha),(3,4,1,2))
	end

	function trade_costs(parameters)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		d = parameters[:d]
		theta = parameters[:theta]
		n_zero = parameters[:numerical_zero]

		kappa = zeros(size(d))
		for j in 1:(J-1)
			for t = 1:T
				kappa[:,:,j,t] = ((d[:,:,j,t] .* transpose(d[:,:,j,t])) ./ (diag(d[:,:,j,t]) * transpose(diag(d[:,:,j,t])))).^(1 / (2 * theta))
			end
		end

		kappa[kappa .< n_zero] = n_zero
		kappa[:,:,end,:] = repeat(eye(N), outer = [1,1,1,T]) # Services

		return kappa = min.(kappa,1)
	end

	function expenditure_shares(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		import_shares = data["import_shares"]
		n_zero = parameters[:numerical_zero]
		
		d = import_shares

		share = d ./ repeat(sum(d,2), outer = [1,N,1,1])
		d_share = 1./repeat(sum(d,2), outer = [1,N,1,1]) - 1
		d_share[d_share .< n_zero] = n_zero
		d = share ./ (1 + d_share)

		for n in 1:N
			d[n,n,:,:] = ones(J,T) - squeeze(sum(d[n,:,:,:],1),1)
		end

		d[d .< n_zero] = n_zero

		return d
	end

	function calculate_xi(parameters)
		theta = parameters[:theta]
		eta = parameters[:eta]

		return gamma((theta + 1 - eta)/theta)
	end

	function calculate_B(parameters)
		beta = parameters[:beta_j]
		gamma = parameters[:gamma_jk]

		gamma = permutedims(cat(ndims(gamma) + 2,gamma),[1,3,2,4])
		return B = (beta .^ -beta) .* prod(gamma .^ -gamma, 1)
	end

	function calculate_psi(parameters, data)
		va = data["va"]
		weights = parameters[:bp_weights]

		va_shares = va ./ repeat(sum(va, 3), outer = [1,1,size(va,3),1])
		_, psi_t = DetrendUtilities.detrend(va_shares, weights)
		return psi_t
	end

	function calculate_p(parameters, data)
		p_sectoral_data = data["p_sectoral_data"]
		d = parameters[:d]
		kappa = parameters[:kappa_mnjt]
		alpha =  parameters[:nu_njt] # Only in case of CD, CES should be different!
		theta = parameters[:theta]

		p_sectoral_base = p_sectoral_data ./ p_sectoral_data[:,:,:,1]

		# US is assumed to be chosen as a base country (US = end), else pwt should be used to do the conversion
		p_sectoral = exp.( mean(1 / theta * log.(d ./ permutedims(cat(ndims(d),d[end,:,:,:]),[4,1,2,3])) - log.(kappa ./ permutedims(cat(ndims(kappa),kappa[end,:,:,:]),[4,1,2,3])), 2) + repeat(permutedims(cat(ndims(p_sectoral_base),log.(p_sectoral_base[:,end,:,:])), [1,4,2,3]), outer = [size(d,1),1,1,1]) )
		temp = Dict{Symbol, Any}()
		temp[:p_sectoral] = p_sectoral

		p_base = permutedims(cat(ndims(p_sectoral_base),prod((p_sectoral_base ./ alpha) .^ alpha, 3)), [1,2,3,4])
		temp[:p_base] = p_base

		p_services = compute_p_services(temp, parameters, data)
		p_sectoral[isnan.(p_sectoral)] = 0
		return p_sectoral = p_sectoral + p_services
	end

	function compute_p_services(temp, parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		p_sectoral = temp[:p_sectoral]
		p_base = temp[:p_base]
		alpha = parameters[:nu_njt] # Only in case of CD, CES should be different!
		pwt = data["pwt"]

		p_services = zeros(N,1,J,T)

		# Service sector is assumed to be on the last position in the sector dimension
		for t in 1:T
			for n in 1:N
				p_services[n,:,end,t] = (pwt[:,n,:,t] * p_base[:,end,:,t]) .^ (1 ./ alpha[:,:,end,t]) .* prod(alpha[:,:,:,t] .^ (-alpha[:,:,:,t])) .^ (-1 ./ alpha[:,:,end,t]) .* prod(p_sectoral[n,:,1:(end-1),t] .^ alpha[1,:,1:(end-1),t]) .^ (-1 ./ alpha[1,:,end,t])
			end
		end

		return p_services
	end

	function calculate_z(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		p_sectoral = parameters[:p_sectoral]
		beta = parameters[:beta_j]
		gamma = parameters[:gamma_jk]
		kappa = parameters[:kappa_mnjt]
		psi = parameters[:psi]
		B = parameters[:B_j]
		d = parameters[:d]
		va = data["va"]
		xi = parameters[:xi]
		theta = parameters[:theta]

		beta = squeeze(beta,(1,2,4))

		z = zeros(1, N, J, T)

		for n in 1:N
			for t in 1:T
				# Sectors except services
				z[1, n, :, t] = exp.( theta * log.(xi * B[1,1,:,1]) + theta * beta .* (log.(va[1,n,:,t]) - log.(psi[1,n,:,t])) + reshape(transpose(mean(log.(d[:,n,:,t]) - theta * log.(kappa[:,n,:,t]),1)) - theta * transpose(mean(log.(p_sectoral[:,1,:,t]),1)),J) + transpose(gamma[:,:,1,1]) * theta * log.(p_sectoral[n,1,:,t]) )

				# Services
				z[1,n,end,t]  = xi^theta * B[1,1,end,1]^theta * (va[1,n,end,t] / psi[1,n,end,t])^(theta * beta[end]) * prod(p_sectoral[n,1,:,t].^(gamma[:,end]))^theta * p_sectoral[n,1,end,t]^(-theta)
			end
		end

		return z
	end

	function calculate_A(parameters)
		z = parameters[:z]
		theta = parameters[:theta]

		return z .^ (1/theta)
	end

	function draw_next_productivity(parameters, t)
		mean_producitivity = exp.(mean(log.(parameters[:A]),4))

		# use "i"th realization to continue future paths
		N, J, S = parameters[:N], parameters[:J], parameters[:S]
		# set variance covariance matrix here
		innovation = exp.(parameters[:shock_stdev] .* randn(1,N,J,S - 1))
		random_realization = non_random_variable(parameters[:A], t)
		AR_decay = parameters[:AR_decay]
		if t > 1
			past_productivity = non_random_variable(parameters[:A], t-1)
			# reversion towards mean
			return cat(4, random_realization, mean_producitivity .^ (1-AR_decay) .* past_productivity .^ (AR_decay) .* innovation)
		else
			# NB: no uncertainty in first period
			return random_realization
		end
	end
end

