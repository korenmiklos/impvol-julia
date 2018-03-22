include("utils.jl")

# for matrix conformity, store all variables in a 4-dimensional array:
# mnjt: destination, source, sector, time

function coerce_parameters!(parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]

	parameters[:alpha_jt] = reshape(parameters[:alpha], (1, 1, J, T))
	parameters[:beta_j] = reshape(parameters[:beta], (1, 1, J, 1))
	parameters[:L_nt] = ones(1,N,1,T)

	beta_j = parameters[:beta]'
	gamma_jk = parameters[:gamma_jk]
	parameters[:B_j_theta] = reshape(beta_j .^ (-beta_j) .* prod(gamma_jk .^ (-gamma_jk),2), (1,1,J,1)) .^ (-parameters[:theta])
end

function coerce_variable!(variable, parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	variable = ones(N,N,J,T) .* variable
end

function price_index(sectoral_prices, globals)
	alpha = globals[:alpha]
	return prod(alpha.^(-alpha) .* sectoral_prices.^(alpha))
end

function input_price_index(sectoral_prices_j, globals)
	gamma_jk = globals[:gamma_jk]
	return prod(sectoral_prices_j .^ gamma_jk, 1)
end

function input_price_index!(variables, parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	variables[:rho_njt] = Array{Float64}(1,N,J,T)
	for n in 1:N
		for t in 1:T
			p = variables[:P_njt][1,n,:,t]
			variables[:rho_njt][1,n,:,t] = input_price_index(p, parameters)
		end
	end
end

function compute_Ds(variables, parameters)
	w_njt = variables[:w_njt]

	Z_njt = parameters[:Z_njt]
	L_nt = parameters[:L_nt]
	beta_j = parameters[:beta_j]
	theta = parameters[:theta]
	kappa_mnjt = parameters[:kappa_mnjt]

	return Z_njt .* (L_nt .* w_njt) .^ (-beta_j*theta) .* kappa_mnjt .^ theta
end

function eq18rhs(variables, parameters)
	input_price_index!(variables, parameters)
	rho_njt = variables[:rho_njt]
	theta = parameters[:theta]
	B_j = parameters[:B_j_theta]
	D = compute_Ds(variables, parameters)
	# FIXME add constant params
	return array_transpose(B_j .* sum(D .* rho_njt .^ (-theta), 2)) .^ (-1/theta)
end

function check_parameters(globals)
	@assert0 sum(globals[:alpha], 2)-1.0
end

parameters = Dict{Symbol, Any}()
N = 24
J = 25
T = 100
fill_dict!(parameters, N=N, J=J, T=T)

alpha = rand(J) .* ones(J,T)
beta = rand(1, J)
kappa_mnjt = rand(N,N,J,T)
for j=1:J
	for t=1:T
		for n=1:N
			kappa_mnjt[n,n,j,t] = 1.0
		end
	end
end
Z_njt = rand(1,N,J,T)

fill_dict!(parameters, alpha=alpha ./ sum(alpha, 1), beta=beta, theta=4, kappa_mnjt=kappa_mnjt, Z_njt=Z_njt)
for (k, v) in parameters
	println(k, size(v))
end
gamma_jk = rand(J,J)
parameters[:gamma_jk] = gamma_jk ./ sum(gamma_jk, 1) .* (1-beta)
coerce_parameters!(parameters)

variables = fill_dict(P_njt=rand(1,N,J,T), w_njt=rand(1,N,J,T))

#check_parameters(parameters)

@time input_price_index!(variables, parameters)
println(size(variables[:rho_njt]))

for k=1:20
	@time new_price = eq18rhs(variables, parameters)
	println(k, ": ", new_price[1,1,1,1])
	variables[:P_njt] = new_price
end


#print(size(compute_Ds(variables, parameters)))