include("utils.jl")
variables = Dict{Symbol, Any}()
parameters = Dict{Symbol, Any}()

# for matrix conformity, store all variables in a 4-dimensional array:
# mnjt: destination, source, sector, time

function coerce_parameters!(parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]

	parameters[:alpha_jt] = reshape(parameters[:alpha], (1, 1, J, T))
	parameters[:beta_j] = reshape(parameters[:beta], (1, 1, J, 1))
	parameters[:L_nt] = ones(1,N,1,T)

	beta_j = parameters[:beta]'
	gamma_jk = parameters[:gamma_jk]
	parameters[:B_j] = reshape(beta_j .^ (-beta_j) .* prod(gamma_jk .^ (-gamma_jk),2), (1,1,J,1)) 
	# FIXME: use definition of xi parameter
	parameters[:xi] = 1 
end

function coerce_variable!(variable, parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	variable = ones(N,N,J,T) .* variable
end

function rotate_sectors(A, y)
	M, N, J, T = size(y)
	K = size(A,1)
	B = zeros(M,N,K,T)
	for m=1:M
		for n=1:N
			for t=1:T
				B[m,n,:,t] = A * y[m,n,:,t] 
			end
		end
	end
	return B
end

function input_price_index(sectoral_prices_j, globals)
	gamma_jk = globals[:gamma_jk]
	return prod(sectoral_prices_j .^ gamma_jk, 1)
end

function input_price_index!(variables, parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	P = variables[:P_njt]
	variables[:input_price_njt] = Array{Float64}(1,N,J,T)
	rho = variables[:input_price_njt]
	for n in 1:N
		for t in 1:T
			p = P[1,n,:,t]
			rho[1,n,:,t] = input_price_index(p, parameters)
		end
	end
end

function compute_price!(variables, parameters)
	theta = parameters[:theta]
	kappa_mnjt = parameters[:kappa_mnjt]
	rho_njt = variables[:rho_njt]

	variables[:P_njt] = array_transpose(sum((rho_njt ./ kappa_mnjt) .^ (-theta), 2) .^ (-1/theta))
end

function compute_price_index!(variables, parameters)
	# so that parameters[:theta] can be referred to as theta
	#@unpack_dictionary variables
	#@unpack_dictionary parameters
	alpha_jt = parameters[:alpha_jt]
	P_njt = variables[:P_njt]

	# use formula on p43 of "paper November 8 2017.pdf"
	variables[:P_nt]=prod(alpha_jt .^ (-alpha_jt) .* P_njt .^ (alpha_jt), 3)
end

function compute_import_shares!(variables, parameters)
	compute_price!(variables, parameters)

	P_mjt = array_transpose(variables[:P_njt])
	theta = parameters[:theta]
	kappa_mnjt = parameters[:kappa_mnjt]
	rho_njt = variables[:rho_njt]

	variables[:d_mnjt] = (rho_njt ./ kappa_mnjt ./ P_mjt) .^ (-theta)
end

function compute_wage!(variables, parameters)
	input_price_index!(variables, parameters)

	P_njt = variables[:P_njt]
	theta = parameters[:theta]
	rho_njt = variables[:rho_njt]
	input_price_njt = variables[:input_price_njt]

	gamma_jk = parameters[:gamma_jk]
	beta_j = parameters[:beta_j]
	B_j = parameters[:B_j]
	xi = parameters[:xi]
	# we dont need T
	A_njt = parameters[:A_njt]

	variables[:w_njt] = (rho_njt .* A_njt ./ input_price_njt ./ (xi * B_j)) .^ (1 ./ beta_j)
end

function compute_revenue!(variables, parameters)
	w_njt = variables[:w_njt]
	L_njt = variables[:L_njt]
	beta_j = parameters[:beta_j]

	variables[:R_njt] = w_njt .* L_njt ./ beta_j
end

function compute_expenditure!(variables, parameters)
	# use eq 19 of "paper November 8 2017.pdf"
	R_nkt = variables[:R_njt]
	beta_j = parameters[:beta_j]
	gamma_jk = parameters[:gamma_jk]
	alpha_jt = parameters[:alpha_jt]
	S_nt = parameters[:S_nt]

	wagebill_nt = rotate_sectors(beta_j[:]',R_nkt)
	intermediate_njt = rotate_sectors(gamma_jk,R_nkt)
	variables[:E_mjt] = array_transpose(alpha_jt .* wagebill_nt + intermediate_njt - alpha_jt .* S_nt)
end

function shadow_price_step(variables, parameters)
	compute_expenditure!(variables, parameters)
	theta = parameters[:theta]
	E_mjt = variables[:E_mjt]
	R_njt = variables[:R_njt]
	kappa_mnjt = parameters[:kappa_mnjt]
	P_mjt = array_transpose(variables[:P_njt])

	return sum( (kappa_mnjt .* P_mjt) .^ theta .* E_mjt ./ R_njt, 1) .^ (1/theta) 
end

function starting_values(variables, parameters)
	compute_expenditure!(variables, parameters)
	theta = parameters[:theta]
	E_mjt = variables[:E_mjt]
	R_njt = variables[:R_njt]
	kappa_mnjt = ones(parameters[:kappa_mnjt])
	P_mjt = array_transpose(variables[:P_njt])

	return sum( (kappa_mnjt .* P_mjt) .^ theta .* E_mjt ./ R_njt, 1) .^ (1/theta) 
end

function loop!(variables, parameters)
	# starting value
	variables[:rho_njt] = starting_values(variables, parameters)

	for k=1:10
		new_rho = shadow_price_step(variables, parameters)
		println(k, mean((new_rho-variables[:rho_njt]) .^ 2))
		variables[:rho_njt] = new_rho
		compute_wage!(variables, parameters)
		compute_price!(variables, parameters)
		compute_revenue!(variables, parameters)
	end
end

function check_parameters(globals)
	@assert0 sum(globals[:alpha], 2)-1.0
end

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

fill_dict!(parameters, alpha=alpha ./ sum(alpha, 1), beta=beta, theta=4, kappa_mnjt=kappa_mnjt, Z_njt=Z_njt, S_nt=zeros(1,N,1,T))
for (k, v) in parameters
	println(k, size(v))
end
gamma_jk = rand(J,J)
parameters[:gamma_jk] = gamma_jk ./ sum(gamma_jk, 1) .* (1-beta)
coerce_parameters!(parameters)

variables = fill_dict(P_njt=rand(1,N,J,T), w_njt=rand(1,N,J,T), R_njt=rand(1,N,J,T))
parameters[:A_njt] = rand(1,N,J,T)
variables[:L_njt] = ones(1,N,J,T)

#check_parameters(parameters)
@time loop!(variables, parameters)


#print(size(compute_Ds(variables, parameters)))