macro assert0(expression)
	:(@assert all($expression .<= 1e-12))
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

