using Images

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

function share_within_year(X)
	return X ./ mean(X, [1 2 3])
end

function distance(p1, p2)
	return meanfinite(((log(p1) .- log(p2)) .^ 2)[:], 1)[1] .^ 0.5
end

function free_trade_sector_shares!(parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	gamma_jk = parameters[:gamma_jk]
	beta_j = parameters[:beta_j]
	alpha_jt = parameters[:alpha_jt]

	revenue_shares = zeros(1,1,J,T)
	for t=1:T
		revenue_shares[1,1,:,t] = eigen_share(alpha_jt[1,1,:,t]*beta_j[1,1,:,1]' + gamma_jk)
	end
	parameters[:sector_shares] = revenue_shares
end

function free_trade_wage_step!(variables, parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	beta_j = parameters[:beta_j]
	theta = parameters[:theta]
	diag_beta = diagm(beta_j[1,1,:,1])
	L_njt = variables[:L_njt]
	A_njt = parameters[:A_njt]
	xi = parameters[:xi]
	B_j = parameters[:B_j]
	P_0t = variables[:P_0t]

	# numeraire: world revenues sum to one in each year t
	# world revenue equals world expenditure
	E_wt = parameters[:revenue_shares]

	RHS = log(beta_j).-log(L_njt).+log(E_wt) .+ theta .* (log(A_njt).-log(xi).-log(B_j).+rotate_sectors(eye(J)-gamma_jk, log(P_0t)))
	ln_w = rotate_sectors(inv(eye(J)+theta*diag_beta), RHS)
	ln_rho = log(xi).+log(B_j).-log(A_njt).+ beta_j.*ln_w .+ rotate_sectors(gamma_jk, log(P_0t))
	new_P = array_transpose(sum(exp(ln_rho).^(-theta), 2) .^ (-1/theta))

	variables[:w_njt] = exp(ln_w)
	variables[:rho_njt] = exp(ln_rho)
	variables[:P_0t] = new_P
end

function free_trade_loop!(variables, parameters)
	compute_free_trade_revenue_shares!(parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	variables[:P_0t] = ones(1,1,J,T) 

	for k=1:20
		previous_price = copy(variables[:P_0t])
		free_trade_wage_step!(variables, parameters)
		println(k, ": ", variables[:w_njt][1,:,:,1])
		println(previous_price[1,1,:,1])
		println(variables[:P_0t][1,1,:,1])
	end
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

function free_trade_country_shares!(variables, parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	A_njt = parameters[:A_njt]
	L_njt = variables[:L_njt]
	theta = parameters[:theta]
	beta_j = parameters[:beta_j]

	d_njt = A_njt .^ (theta ./ (1 + theta*beta_j)) .* L_njt .^ (theta .* beta_j ./ (1 + theta*beta_j))
	variables[:d_njt_free] = d_njt ./ sum(d_njt,2)
end

function free_trade_wages!(variables, parameters)
	free_trade_country_shares!(variables, parameters)
	free_trade_sector_shares!(parameters)

	d_njt = variables[:d_njt_free]
	beta_j = parameters[:beta_j]
	E_wt = parameters[:sector_shares]
	L_njt = variables[:L_njt]

	variables[:w_njt] = beta_j .* d_njt .* E_wt ./ L_njt
end

function free_trade_prices!(variables, parameters)
	free_trade_country_shares!(variables, parameters)
	free_trade_sector_shares!(parameters)
	beta_j = parameters[:beta_j]
	theta = parameters[:theta]
	d_njt = variables[:d_njt_free]
	xi = parameters[:xi]
	B_j = parameters[:B_j]
	A_njt = parameters[:A_njt]
	E_wt = parameters[:sector_shares]
	L_njt = variables[:L_njt]
	gamma_jk = parameters[:gamma_jk]

	variables[:P_njt] = exp(rotate_sectors(inv(eye(gamma_jk)-gamma_jk), log(xi * d_njt .^(beta_j+1/theta) .* B_j .* (E_wt ./ L_njt) .^(beta_j) ./ A_njt)))
	variables[:rho_njt] = variables[:P_njt] ./ d_njt .^(1/theta)
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

	rho_njt = variables[:rho_njt]
	input_price_njt = variables[:input_price_njt]

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
	theta = parameters[:theta]
	E_mjt = variables[:E_mjt]
	R_njt = variables[:R_njt]
	kappa_mnjt = parameters[:kappa_mnjt]
	P_mjt = array_transpose(variables[:P_njt])

	return sum( (kappa_mnjt .* P_mjt) .^ theta .* E_mjt ./ R_njt, 1) .^ (1/theta)  
end

function loop!(variables, parameters)
	# starting value
	lambda = parameters[:lambda]
	dist = 999
	k = 1

	while dist > parameters[:tolerance]
		new_rho = shadow_price_step(variables, parameters)
		dist = distance(new_rho, variables[:rho_njt])
		println(k, ": ")
		println(meanfinite(new_rho, 4)[1,1,1,1])
		println(dist)
		variables[:rho_njt] = lambda*new_rho + (1-lambda)*variables[:rho_njt]
		compute_price!(variables, parameters)
		compute_wage!(variables, parameters)
		compute_revenue!(variables, parameters)
		compute_expenditure!(variables, parameters)
		k = k+1
	end
end

function check_parameters(globals)
	@assert0 sum(globals[:alpha], 2)-1.0
end

N = 24
J = 25
T = 1000
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
gamma_jk = rand(J,J)
# for testing purposed, set IO links to 0
parameters[:gamma_jk] = zeros(J,J)
parameters[:beta] = ones(1,J)
#parameters[:gamma_jk] = gamma_jk ./ sum(gamma_jk, 1) .* (1-beta)
# adaptive step size. large lambda means large steps
parameters[:lambda] = exp(-0.10*(J-1)^0.75)
# this is log points of average input price differences
parameters[:tolerance] = 0.001
coerce_parameters!(parameters)

variables = fill_dict(P_njt=rand(1,N,J,T), w_njt=rand(1,N,J,T), R_njt=rand(1,N,J,T))
parameters[:A_njt] = rand(1,N,J,T)
variables[:L_njt] = ones(1,N,J,T)

free_trade_wages!(variables, parameters)
free_trade_prices!(variables, parameters)
compute_revenue!(variables, parameters)
compute_expenditure!(variables, parameters)


println(variables[:d_njt_free][1,:,1,1])
println(variables[:w_njt][1,:,1,1])
println(variables[:P_njt][1,:,:,1])
println(variables[:rho_njt][1,:,:,1])

@time loop!(variables, parameters)

