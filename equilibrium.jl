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
	parameters[:B_j_theta] = reshape(beta_j .^ (-beta_j) .* prod(gamma_jk .^ (-gamma_jk),2), (1,1,J,1)) .^ (-parameters[:theta])
end

function coerce_variable!(variable, parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	variable = ones(N,N,J,T) .* variable
end

function rotate_sectors(A, y)
	N, M, J, T = size(y)
	K = size(A,1)
	B = zeros(N,M,K,T)
	for n=1:N
		for m=1:M
			for t=1:T
				B[n,m,:,t] = A * y[n,m,:,t] 
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
	variables[:rho_njt] = Array{Float64}(1,N,J,T)
	rho = variables[:rho_njt]
	for n in 1:N
		for t in 1:T
			p = P[1,n,:,t]
			rho[1,n,:,t] = input_price_index(p, parameters)
		end
	end
end

function compute_Ds!(variables, parameters)
	w_njt = variables[:w_njt]
	Z_njt = variables[:Z_njt]
	L_nt = parameters[:L_nt]
	beta_j = parameters[:beta_j]
	theta = parameters[:theta]
	kappa_mnjt = parameters[:kappa_mnjt]
	# use formula on p43 of "paper November 8 2017.pdf"
	variables[:D] = Z_njt .* (L_nt .* w_njt) .^ (-beta_j*theta) .* kappa_mnjt .^ theta
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
	compute_Ds!(variables, parameters)
	compute_price_index!(variables, parameters)

	D = variables[:D]
	P_nt = variables[:P_nt]
	theta = parameters[:theta]
	beta_j = parameters[:beta_j]

	# use formula on p45 of "paper November 8 2017.pdf"
	share = D .* (P_nt .^ (theta*(beta_j-1)))
	variables[:d_mnjt] = share ./ sum(share,2)
end

function price_step(variables, parameters)
	input_price_index!(variables, parameters, true)

	B_j_theta = parameters[:B_j_theta]
	D = variables[:D]
	rho_njt = variables[:rho_njt_theta]
	theta = parameters[:theta]

	# FIXME add constant params
	return array_transpose(B_j_theta .* sum(D .* rho_njt .^ (-theta), 2)) .^ (-1/theta)
end

function price_loop!(variables, parameters)
	compute_Ds!(variables, parameters)
	# FIXME: init with a reasobale price vector
	for k=1:30
		# FIXME: check convergence
		new_price = price_step(variables, parameters)
		variables[:P_njt_theta] = new_price
	end
end

function revenue_step(variables, parameters)
	compute_import_shares!(variables, parameters)

	d_mnjt = variables[:d_mnjt]
	alpha_jt = parameters[:alpha_jt]
	beta_j = parameters[:beta_j]
	R_njt = variables[:R_njt]
	S_nt = parameters[:S_nt]

	# use eq 19 of "paper November 8 2017.pdf"
	wagebill_nt = rotate_sectors(beta_j[:]',R_njt)
	intermediate_njt = rotate_sectors(gamma_jk,R_njt)
	return sum(array_transpose(d_mnjt) .* (alpha_jt .* wagebill_nt + intermediate_njt - alpha_jt .* S_nt), 1)
end

function revenue_loop!(variables, parameters)
	# FIXME: init with a reasobale revenue vector
	for k=1:30
		# FIXME: check convergence
		new_revenue = revenue_step(variables, parameters)
		println(new_revenue[1,1])
		variables[:R_njt] = new_revenue
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

global variables = fill_dict(P_njt=rand(1,N,J,T), w_njt=rand(1,N,J,T), Z_njt=rand(1,N,J,T), R_njt=rand(1,N,J,T))

#check_parameters(parameters)
@time input_price_index!(variables, parameters)
println(size(variables[:rho_njt]))

@time price_loop!(variables, parameters)

@time revenue_loop!(variables, parameters)


#print(size(compute_Ds(variables, parameters)))