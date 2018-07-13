function manipulate_import_shares(data::DataFrames.DataFrame,dim_size::NTuple)
	# This function extends the imported import shares, so that it will be full dimensional
	# Adds service sector (with 0 import share) and adds own country import share (0)
	insert!(data, length(names(data))  + 1, 0, :s)
	N = dim_size[2]
	J = dim_size[4] + 1
	T = dim_size[3]
	for t in 1:T
		y = 1971 + t
		for d in 1:N
			vec = zeros(1,J + 3)
			vec[1,1:3] = [y,d,d]
			push!(data, vec)
		end
	end
	return sort!(data,(1,2,3)), map(+,dim_size,(1,0,0,1))
end

function compute_gammas(beta::Array{Float64,4},io_values::Array{Float64,4},total_output::Array{Float64,4},output_shares::Array{Float64,4},intermediate_input_shares::Array{Float64,4})
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
	output_shares_full[split_idx,:,:,:] = permutedims(repeat(output_shares,outer = [1,1,size(io_values_dupl,2),1]),(2,3,1,4))

	io_values_new = io_values_dupl .* output_shares_full

	# Split columns
	io_values_dupl = io_values_new[:,dupl_idx,:,:]

	intermediate_input_shares_full = ones(size(io_values_dupl))
	intermediate_input_shares_full[:,split_idx,:,:] = repeat(intermediate_input_shares,outer = [size(io_values_dupl,1),1,1,1])

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
	gammas = io_values_new ./ repeat(total_output,outer = [size(io_values_new,1),1,1,1])
	gammas = mean(gammas,4)
	gammas = gammas .* permutedims(1-beta,(1,3,2,4)) ./ sum(gammas,1)
	return gammas = repeat(gammas, outer = [1,1,1,T])
end

# function compute_gammas(beta::Array{Float64,2},io_values::Array{Float64,3},total_output::Array{Float64,2},output_shares::Array{Float64,2},intermediate_input_shares::Array{Float64,2})
# 	# Summing sectors
# 	# Agriculture & mining
# 	io_values[1,:,:]  = sum(io_values[1:2,:,:],1)
# 	io_values         = io_values[setdiff(1:end,2),:,:]
# 	io_values[:,1,:]  = sum(io_values[:,1:2,:],2)
# 	io_values         = io_values[:,setdiff(1:end,2),:]
# 	total_output[1,:] = sum(total_output[1:2,:],1)
# 	total_output      = total_output[setdiff(1:end,2),:]

# 	# Services
# 	io_values[18,:,:]  = sum(io_values[18:end,:,:],1)
# 	io_values          = io_values[setdiff(1:end,19:end),:,:]
# 	io_values[:,18,:]  = sum(io_values[:,18:end,:],2)
# 	io_values          = io_values[:,setdiff(1:end,19:end),:]
# 	total_output[18,:] = sum(total_output[18:end,:],1)
# 	total_output       = total_output[setdiff(1:end,19:end),:]

# 	# Split rows (I don't know what's going on)
# 	dupl_idx = [1, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 14, 15, 16, 17, 18]
# 	io_values_dupl = io_values[dupl_idx,:,:]

# 	output_shares_full = ones(size(io_values_dupl))

# 	split_idx = [2, 3, 4, 5, 6, 8, 9, 17, 18, 19]
# 	output_shares_full[split_idx,:,:] = permutedims(repeat(output_shares,outer = [1,1,size(io_values_dupl,2)]),(1,3,2))

# 	io_values_new = io_values_dupl .* output_shares_full

# 	# Split columns (I still don't know what's going on)
# 	io_values_dupl = io_values_new[:,dupl_idx,:]

# 	intermediate_input_shares_full = ones(size(io_values_dupl))
# 	intermediate_input_shares_full[:,split_idx,:] = permutedims(repeat(intermediate_input_shares,outer = [1,1,size(io_values_dupl,1)]),(3,1,2))

# 	io_values_new = io_values_dupl .* intermediate_input_shares_full

# 	total_output = total_output[dupl_idx,:]
# 	output_shares_full = ones(size(total_output))
# 	output_shares_full[split_idx,:] = output_shares
# 	total_output = total_output .* output_shares_full

# 	# Correct order of sectors
# 	io_values_new[18,:,:], io_values_new[19,:,:], io_values_new[20,:,:] = io_values_new[20,:,:], io_values_new[18,:,:], io_values_new[19,:,:]
# 	io_values_new[:,18,:], io_values_new[:,19,:], io_values_new[:,20,:] = io_values_new[:,20,:], io_values_new[:,18,:], io_values_new[:,19,:]
# 	total_output[18,:], total_output[19,:], total_output[20,:] = total_output[20,:], total_output[18,:], total_output[19,:]

# 	gammas = io_values_new ./ permutedims(repeat(total_output,outer = [1,1,size(io_values_new,1)]),(3,1,2))
# 	gammas = squeeze(mean(gammas,3),3)
# 	return gammas = gammas .* array_transpose(1-beta) ./ sum(gammas,1)
# end

function compute_alphas(va::Array{Float64,4}, beta::Array{Float64,4}, gammas::Array{Float64,4},weights::Array{Float64,1})
	M, N, J, T = size(beta)

	beta = mean(beta,(1,2,4))

	alphas = zeros(J,T)

	va = squeeze(va,1)
	beta = squeeze(beta,(1,4))
	gammas = squeeze(mean(gammas,4),(3,4))

	for t in 1:T
		va_t = transpose(sum(va[:,:,t],1))
		alphas[:,t] = (eye(J) - gammas) * diagm(1 ./ beta[:],0) * va_t / sum(va_t)
	end

	# Replace negative elements with 0
	alphas = (alphas + abs.(alphas)) / 2

	# Smooth the series
	alpha_c, alpha_t = detrend(alphas, weights)

	# Normalization
	alphas = alpha_t ./ sum(alpha_t,1)

	return alphas = permutedims(cat(ndims(alphas) + 2,alphas),(3,4,1,2))
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

function trade_costs(d::Array{Float64,4}, theta::Float64, n_zero::Float64)
	M, N, J, T = size(d)

	kappa = zeros(size(d))
	for j in 1:(J-1)
		for t = 1:T
			kappa[:,:,j,t] = ((d[:,:,j,t] .* transpose(d[:,:,j,t])) ./ (diag(d[:,:,j,t]) * transpose(diag(d[:,:,j,t])))).^(1 / (2 * theta))
		end
	end

	kappa[:,:,end,:] = repeat(eye(M), outer = [1,1,1,T])

	kappa[kappa .< n_zero] = n_zero
	kappa = min.(kappa,1)
end