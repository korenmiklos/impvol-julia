
module CalibrateParameters
	include("calibration_utils.jl")
	using JLD
	
	function calibrate_parameters!(parameters)
		data = JLD.load("../../../data/impvol_data.jld")

		# standard deviation for each (n,j)
		parameters[:shock_stdev] = 0.1*ones(1,N,J,1)
		# AR coefficient for each (n,j)
		parameters[:AR_decay] = 0.9*ones(1,N,J,1)

		parameters[:beta] = data["beta"]
		_, N, J, T = size(parameters[:beta])
		parameters[:N], parameters[:J], parameters[:T] = N, J, T

		parameters[:gamma_jk] = compute_gammas(parameters[:beta],data["io_values"],data["total_output"],data["output_shares"],data["intermediate_input_shares"])
		# CD case
		parameters[:nu_njt] = compute_alphas(data["va"],parameters[:beta],parameters[:gamma_jk],parameters[:bp_weights])
 
		parameters[:S_nt] = zeros(1,N,1,T)

		d = expenditure_shares(data["import_shares"], parameters[:numerical_zero])

		parameters[:kappa] = trade_costs(d, parameters[:theta], parameters[:numerical_zero])

		p_sectoral = calculate_p(data["p_sectoral_data"], data["pwt"], d, parameters[:kappa], parameters[:nu_njt], parameters[:theta])
		psi = calculate_psi(data["va"], parameters[:bp_weights])
		B = calculate_B(parameters[:beta], parameters[:gamma_jk])
		xi = calculate_xi(parameters[:theta], parameters[:eta])
		# FIXME: pls save these into parameters[]

		parameters[:z] = calculate_z(p_sectoral, parameters[:beta], parameters[:gamma_jk], parameters[:kappa], psi, B, d, data["va"], xi, parameters[:theta])
	end

	function compute_gammas(beta::Array{Float64,4}, io_values::Array{Float64,4}, total_output::Array{Float64,4}, output_shares::Array{Float64,4}, intermediate_input_shares::Array{Float64,4})
		M, N, J, T = size(beta)
		beta = mean(beta,(1,2,4))

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

		# Compute gammas
		gammas = io_values_new ./ repeat(total_output, outer = [size(io_values_new,1),1,1,1])
		gammas = mean(gammas,4)
		gammas = gammas .* permutedims(1-beta,(1,3,2,4)) ./ sum(gammas,1)
		return gammas = squeeze(gammas,(3,4))
	end

	function compute_alphas(va::Array{Float64,4}, beta::Array{Float64,4}, gammas::Array{Float64,2},weights::Array{Float64,1})
		M, N, J, T = size(beta)

		beta = mean(beta,(1,2,4))

		alphas = zeros(J,T)

		va = squeeze(va,1)
		beta = squeeze(beta, (1,4))

		for t in 1:T
			va_t = transpose(sum(va[:,:,t],1))
			alphas[:,t] = (eye(J) - gammas) * diagm(1 ./ beta[:],0) * va_t / sum(va_t)
		end

		# Replace negative elements with 0
		alphas = (alphas + abs.(alphas)) / 2

		# Smooth the series
		alpha_c, alpha_t = DetrendUtilities.detrend(alphas, weights)

		# Normalization
		alphas = alpha_t ./ sum(alpha_t,1)

		return alphas = permutedims(cat(ndims(alphas) + 2,alphas),(3,4,1,2))
	end

	function trade_costs(d::Array{Float64,4}, theta::Float64, n_zero::Float64)
		M, N, J, T = size(d)

		kappa = zeros(size(d))
		for j in 1:(J-1)
			for t = 1:T
				kappa[:,:,j,t] = ((d[:,:,j,t] .* transpose(d[:,:,j,t])) ./ (diag(d[:,:,j,t]) * transpose(diag(d[:,:,j,t])))).^(1 / (2 * theta))
			end
		end

		kappa[kappa .< n_zero] = n_zero
		kappa[:,:,end,:] = repeat(eye(M), outer = [1,1,1,T]) # Services

		return kappa = min.(kappa,1)
	end

	function expenditure_shares(import_shares::Array{Float64,4}, n_zero::Float64)
		M, N, J, T = size(import_shares)
		d = import_shares

		share = d ./ repeat(sum(d,2), outer = [1,N,1,1])
		d_share = 1./repeat(sum(d,2), outer = [1,N,1,1]) - 1
		d_share[d_share .< n_zero] = n_zero
		d = share ./ (1 + d_share)

		for m in 1:M
			d[m,m,:,:] = ones(J,T) - squeeze(sum(d[m,:,:,:],1),1)
		end

		d[d .< n_zero] = n_zero

		return d
	end

	function calculate_xi(theta::Float64, eta::Float64)
		return gamma((theta + 1 - eta)/theta)
	end

	function calculate_B(beta::Array{Float64,4}, gammas::Array{Float64,2})
		gammas = cat(ndims(gammas) + 2,gammas)
		beta = squeeze(mean(beta,(1,2,4)),(1,4))
		return B = (beta .^ -beta) .* permutedims(prod(gammas .^ -gammas, 1), [2,1])
	end

	function calculate_psi(va::Array{Float64,4}, weights::Array{Float64,1})
		va_shares = va ./ repeat(sum(va, 3), outer = [1,1,size(va,3),1])
		_, psi_t = DetrendUtilities.detrend(va_shares, weights)
		return psi_t
	end

	function calculate_p(p_sectoral_data::Array{Float64,4}, pwt::Array{Float64,4}, d::Array{Float64,4}, kappa::Array{Float64,4}, alphas::Array{Float64,4}, theta::Float64)
		p_sectoral_base = p_sectoral_data ./ p_sectoral_data[:,:,:,1]

		# US is assumed to be chosen as a base country (US = end), else pwt should be used to do the conversion
		p_sectoral = exp.( mean(1 / theta * log.(d ./ permutedims(cat(ndims(d),d[end,:,:,:]),[4,1,2,3])) - log.(kappa ./ permutedims(cat(ndims(kappa),kappa[end,:,:,:]),[4,1,2,3])), 2) + repeat(permutedims(cat(ndims(p_sectoral_base),log.(p_sectoral_base[:,end,:,:])), [1,4,2,3]), outer = [size(d,1),1,1,1]) )

		p_base = permutedims(cat(ndims(p_sectoral_base),prod((p_sectoral_base ./ alphas) .^ alphas, 3)), [1,2,3,4])
		p_services = compute_p_services(p_sectoral, p_base, pwt, alphas)
		p_sectoral[isnan.(p_sectoral)] = 0
		return p_sectoral = p_sectoral + p_services
	end

	function compute_p_services(p_sectoral::Array{Float64,4}, p_base::Array{Float64,4}, pwt::Array{Float64,4}, alphas::Array{Float64,4})
		M, N, J, T = size(p_sectoral)
		p_services = zeros(M,N,J,T)

		# Service sector is assumed to be on the last position in the sector dimension
		for t in 1:T
			for m in 1:M
				p_services[m,:,end,t] = (pwt[:,m,:,t] * p_base[:,end,:,t]) .^ (1 ./ alphas[:,:,end,t]) .* prod(alphas[:,:,:,t] .^ (-alphas[:,:,:,t])) .^ (-1 ./ alphas[:,:,end,t]) .* prod(p_sectoral[m,:,1:(end-1),t] .^ alphas[1,:,1:(end-1),t]) .^ (-1 ./ alphas[1,:,end,t])
			end
		end

		return p_services
	end

	function calculate_z(p_sectoral::Array{Float64,4}, beta::Array{Float64,4}, gammas::Array{Float64,2}, kappa::Array{Float64,4}, psi::Array{Float64,4}, B::Array{Float64,2}, d::Array{Float64,4}, va::Array{Float64,4}, xi::Float64, theta::Float64)
		_, N, J, T = size(va)
		beta = squeeze(mean(beta,(1,2,4)),(1,2,4))

		z = zeros(1, N, J, T)

		for n in 1:N
			for t in 1:T
				# Sectors except services
				z[1, n, :, t] = exp.( theta * log.(xi * B[1,1,:,1]) + theta * beta .* (log.(va[1,n,:,t]) - log.(psi[1,n,:,t])) + reshape(transpose(mean(log.(d[:,n,:,t]) - theta * log.(kappa[:,n,:,t]),1)) - theta * transpose(mean(log.(p_sectoral[:,1,:,t]),1)),J) + transpose(gammas[:,:,1,1]) * theta * log.(p_sectoral[n,1,:,t]) )

				# Services
				z[1,n,end,t]  = xi^theta * B[1,end]^theta * (va[1,n,end,t] / psi[1,n,end,t])^(theta * beta[end]) * prod(p_sectoral[n,1,:,t].^(gammas[:,end]))^theta * p_sectoral[n,1,end,t]^(-theta)
			end
		end

		return z
	end

	function calculate_A(z::Array{Float64,4}, theta::Float64)
		return z .^ (1/theta)
	end
end

