macro assert0(expression)
	:(@assert all($expression .<= 1e-12))
end

macro unpack_dictionary(expression)
	lines = []
	dict = eval(expression)
	for (key, value) in dict
		:(local $key = $value)
	end 
end

function fill_dict(;kwargs...)
	return Dict(kwargs)
end

function fill_dict!(dict; kwargs...)
	merge!(dict, Dict(kwargs))
end

#SUM_j = kron(eye(N), ones(J,J))
#SUM_n = kron(ones(N,N), eye(J)) 
function sum_j(variable)
	# sum over sectors in a stacked vector/matrix
#	return SUM_j*variable
end
function sum_n(variable)
	# sum over countries in a stacked vector/matrix
#	return SUM_n*variable
end
function prod_j(varialble)
	
end

function array_transpose(array)
	index = [i for i in 1:ndims(array)]
	index[1], index[2] = 2, 1
	return permutedims(array, index)	
end

function rotate_along_dimension(rotator, array, n)
	index = [i for i in 1:ndims(array)]
	index[1], index[n] = n, 1
	temp = permutedims(array, index)
	S1 = size(temp)
	A = zeros()
	for i in CartesianRange(S1[2:end])
		@show i
		A[:,i] = rotator*temp[:,i]
	end
	return A
end

function share_within_dimension(X, i)
	return X ./ sum(X,i)
end

function eigen_share(A)
	# the solution to x = Ax with x'1 = 1
	D, V = eig(A)
	i = find(abs.(D - 1.0) .< 1e-6)[1]
	return V[:,i] ./ sum(V[:,i])
end

function non_random_variable(y, t)
	# array coersion is going to take care of rest
	B = y[:,:,:,t]
	return cat(ndims(B)+1, B)
end
