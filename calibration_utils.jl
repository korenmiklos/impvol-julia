function detrend(X::Array{Float64,2}, weights::Array{Float64,1})
	N, T = size(X)

	X_c = zeros(N,T)
	X_t = zeros(N,T)

	# Maximum lag
	K = size(weights,1) - 1
	w = [flipdim(weights,1); weights[2:end]]

	# Extended series with K 0's at both ends
	X_ext = zeros(N, T + 2K)
	X_ext[:,(K + 1):(K + T)] = X

	# Reflect head/tail
	X_ext[:,1:K] = repeat(X[:,1], outer = [1,K]) + flipdim(repeat(X[:,1], outer = [1,K]) - X[:,2:(K + 1)],2)
	X_ext[:,(T + K + 1):(T + 2K)] = repeat(X[:,T], outer = [1,K]) + flipdim(repeat(X[:,T],outer = [1,K]) - X[:,(T - K):(T - 1)],2)

	for t in 1:T
		X_c[:,t] = X_ext[:,t:(t + 2K)] * w
		X_t[:,t] = X_ext[:,t + K] - X_c[:,t]
	end

	return X_c, X_t
end