module ImpvolEquilibrium

export period_wrapper, coerce_parameters!

using Logging

include("utils.jl")
parameters = Dict{Symbol, Any}()

# for matrix conformity, store all variables in a 4-dimensional array:
# mnjt: destination, source, sector, time

# per-period random variables are stored as
# mnjs: destination, source, sector, state

function expected_value(y)
	return mean(y, 4)[:,:,:,1]
end

function rotate_sectors(A, y)
	# y may be a time series variable or a random variable
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

function distance(p1, p2)
	# distance should only depend on real prices, not on nominal
	normalized_p1 = p1 ./ mean(p1, (1, 2, 3))
	normalized_p2 = p2 ./ mean(p2, (1, 2, 3))
	return mean(((log.(normalized_p1) .- log.(normalized_p2)) .^ 2)[:], 1)[1] .^ 0.5
end

function free_trade_sector_shares!(parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	gamma_jk = parameters[:gamma_jk]
	beta_j = parameters[:beta_j]
	alpha_jt = parameters[:nu_njt] ./ sum(parameters[:nu_njt], 3)

	revenue_shares = zeros(1,1,J,T)
	for t=1:T
		revenue_shares[1,1,:,t] = eigen_share(alpha_jt[1,1,:,t]*beta_j[1,1,:,1]' + gamma_jk)
	end
	# FIXME: choose units such that world expenditure equals that in data
	parameters[:sector_shares] = revenue_shares
end

function input_price_index(sectoral_prices_j, globals)
	gamma_jk = globals[:gamma_jk]
	return prod(sectoral_prices_j .^ gamma_jk, 1)
end

function input_price_index!(random_variables, parameters)
	P = random_variables[:P_njs]
	_, N, J, S = size(P)
	random_variables[:input_price_njs] = Array{Float64}(1,N,J,S)
	rho = random_variables[:input_price_njs]
	for n in 1:N
		for s in 1:S
			p = P[1,n,:,s]
			rho[1,n,:,s] = input_price_index(p, parameters)
		end
	end
end

function compute_price!(random_variables, parameters, t)
	theta = parameters[:theta]
	kappa = non_random_variable(parameters[:kappa_mnjt], t)
	rho_njs = random_variables[:rho_njs]

	random_variables[:P_njs] = array_transpose(sum((rho_njs ./ kappa) .^ (-theta), 2) .^ (-1/theta))
end

function compute_price_index!(random_variables, parameters, t)
	nu = non_random_variable(parameters[:nu_njt], t)
	alpha = nu ./ sum(nu, 3)
	P_njs = random_variables[:P_njs]
	sigma = parameters[:sigma]

	# use formula on p43 of "paper November 8 2017.pdf"
	# Cobb-Douglas is a special case when sigma ~ 1
	random_variables[:P_ns] = sum(alpha .* P_njs .^ (1-sigma), 3) .^ (1/(1-sigma))
end

function free_trade_country_shares!(random_variables, parameters)
	A_njs = random_variables[:A_njs]
	L_njs = random_variables[:L_njs]
	theta = parameters[:theta]
	beta_j = parameters[:beta_j]

	d_njs = A_njs .^ (theta ./ (1 + theta*beta_j)) .* L_njs .^ (theta .* beta_j ./ (1 + theta*beta_j))
	random_variables[:d_njs_free] = d_njs ./ sum(d_njs,2)
end

function free_trade_wages!(random_variables, parameters, t)
	free_trade_country_shares!(random_variables, parameters)
	free_trade_sector_shares!(parameters)

	d_njs = random_variables[:d_njs_free]
	L_njs = random_variables[:L_njs]
	beta_j = parameters[:beta_j]
	E_wt = non_random_variable(parameters[:sector_shares], t) 

	random_variables[:w_njs] = beta_j .* d_njs .* E_wt ./ L_njs
end

function free_trade_prices!(random_variables, parameters, t)
	free_trade_country_shares!(random_variables, parameters)
	free_trade_sector_shares!(parameters)
	beta_j = parameters[:beta_j]
	theta = parameters[:theta]
	d_njs = random_variables[:d_njs_free]
	L_njs = random_variables[:L_njs]
	xi = parameters[:xi]
	B_j = parameters[:B_j]
	A_njs = random_variables[:A_njs]
	E_wt = non_random_variable(parameters[:sector_shares], t)
	gamma_jk = parameters[:gamma_jk]

	random_variables[:P_njs] = exp.(rotate_sectors(inv(eye(gamma_jk)-gamma_jk), log.(xi * d_njs .^(beta_j+1/theta) .* B_j .* (E_wt ./ L_njs) .^(beta_j) ./ A_njs)))
	random_variables[:rho_njs] = random_variables[:P_njs] ./ d_njs .^(1/theta)
end

function compute_wage!(random_variables, parameters)
	input_price_index!(random_variables, parameters)

	rho_njs = random_variables[:rho_njs]
	input_price_njs = random_variables[:input_price_njs]

	beta_j = parameters[:beta_j]
	B_j = parameters[:B_j]
	xi = parameters[:xi]
	# we dont need T
	A_njs = random_variables[:A_njs]

	random_variables[:w_njs] = (rho_njs .* A_njs ./ input_price_njs ./ (xi * B_j)) .^ (1 ./ beta_j)
end

function compute_revenue!(random_variables, parameters)
	w_njs = random_variables[:w_njs]
	L_njs = random_variables[:L_njs]
	beta_j = parameters[:beta_j]

	random_variables[:R_njs] = w_njs .* L_njs ./ beta_j
end

function CES_share(nu, price, sigma)
	temp = nu .* price .^ (1-sigma)
	return temp ./ sum(temp, 3)
end

function compute_expenditure_shares!(random_variables, parameters, t)
	# use eq 19 of "paper November 8 2017.pdf"
	R_nks = random_variables[:R_njs]
	beta_j = parameters[:beta_j]
	gamma_jk = parameters[:gamma_jk]
	nu = non_random_variable(parameters[:nu_njt], t)
	# encompass CES and Cobb-Douglas
	alpha_njt = CES_share(nu, random_variables[:P_njs], parameters[:sigma])
	S_nt = non_random_variable(parameters[:S_nt], t)
	expenditure = sum(R_nks, 3) .- S_nt

	wagebill_ns = rotate_sectors(beta_j[:]', R_nks)
	intermediate_njs = rotate_sectors(gamma_jk, R_nks)
	random_variables[:e_mjs] = array_transpose((alpha_njt .* wagebill_ns .+ intermediate_njs .- alpha_njt .* S_nt) ./ expenditure)
	debug(alpha_njt[1,1,1:2,1])
end

function compute_real_gdp!(random_variables, parameters, t)
	compute_price_index!(random_variables, parameters, t)
	w_njs = random_variables[:w_njs]
	L_njs = random_variables[:L_njs]
	P_ns = random_variables[:P_ns]

	random_variables[:real_GDP] = w_njs .* L_njs ./ P_ns
end

function fixed_expenditure_shares!(random_variables, parameters, t)
	R_mjs = array_transpose(random_variables[:R_njs])	
	expenditure = sum(R_mjs, 3) .- array_transpose(non_random_variable(parameters[:S_nt], t))
	random_variables[:E_mjs] = random_variables[:e_mjs] .* expenditure
end

function shadow_price_step(random_variables, parameters, t)
	theta = parameters[:theta]
	E_mjs = random_variables[:E_mjs]
	R_njs = random_variables[:R_njs]
	kappa_mnjt = non_random_variable(parameters[:kappa_mnjt], t)
	P_mjs = array_transpose(random_variables[:P_njs])

	return sum( (kappa_mnjt .* P_mjs) .^ theta .* E_mjs ./ R_njs, 1) .^ (1/theta)  
end

function starting_values!(random_variables, parameters, t)
	free_trade_wages!(random_variables, parameters, t)
	free_trade_prices!(random_variables, parameters, t)
	compute_revenue!(random_variables, parameters)
	compute_expenditure_shares!(random_variables, parameters, t)
	fixed_expenditure_shares!(random_variables, parameters, t)
end

function inner_loop!(random_variables, parameters, t)
	debug("------ BEGIN Inner loop")
	lambda = parameters[:inner_step_size]
	dist = 999
	k = 1

	while (dist > parameters[:inner_tolerance]) && (k <= parameters[:max_iter_inner])
		new_rho = shadow_price_step(random_variables, parameters, t)
		dist = distance(new_rho, random_variables[:rho_njs])
		debug("-------- Inner ", k, ": ", dist)
		random_variables[:rho_njs] = lambda*new_rho + (1-lambda)*random_variables[:rho_njs]
		compute_price!(random_variables, parameters, t)
		compute_wage!(random_variables, parameters)
		compute_revenue!(random_variables, parameters)
		fixed_expenditure_shares!(random_variables, parameters, t)
		k = k+1
	end
	#warn("inner: ", k-1)
	debug("------ END Inner loop")
end

function middle_loop!(random_variables, parameters, t)
	debug("---- BEGIN Middle loop")
	dist = 999
	k = 1
	old_expenditure_shares = random_variables[:e_mjs]

	while (dist > parameters[:middle_tolerance]) && (k <= parameters[:max_iter_middle])
		inner_loop!(random_variables, parameters, t)
		compute_expenditure_shares!(random_variables, parameters, t)
		dist = distance(random_variables[:e_mjs], old_expenditure_shares)

		random_variables[:e_mjs] = parameters[:middle_step_size]*random_variables[:e_mjs]+(1-parameters[:middle_step_size])*old_expenditure_shares
		info("------ Middle ", k, ": ", dist)

		old_expenditure_shares = random_variables[:e_mjs]
		k = k+1
	end
	#warn("middle: ", k-1)
	debug("---- END Middle loop")
end

function adjustment_loop!(random_variables, L_nj_star, parameters, t)

	debug("-- BEGIN Adjustment loop")
	L_njs = L_nj_star
	L_n = sum(L_nj_star, 3)
	dist = 999
	k = 1
	step_size = parameters[:adjustment_step_size]
	# FIXME: DEBUG: no adjustment
	#step_size = 0.0

	nulla = ones(1,1,1,1)*parameters[:numerical_zero]

	while (dist > parameters[:adjustment_tolerance]) && (k <= parameters[:max_iter_adjustment])
		random_variables[:L_njs] = L_njs
		middle_loop!(random_variables, parameters, t)
		w_njs = random_variables[:w_njs]
		w_ns = sum(w_njs .* L_njs, 3) ./ sum(L_njs, 3)
		wage_gap = w_njs ./ w_ns

		# BUGFIX: lambda L_nt is not always 1, lambda is such that epsilon sums to zero
		epsilon_per_L = parameters[:one_over_rho]*wage_gap
		epsilon_per_L = epsilon_per_L .- mean(epsilon_per_L, 3)
		debug("Adjustment: ", epsilon_per_L[1,1,:,1])
		epsilon_per_L = max.(nulla, epsilon_per_L)
		epsilon_per_L = min.(1 .- nulla, epsilon_per_L)
		debug("Min wage gap: ", minimum(wage_gap))
		debug("Max wage gap: ", maximum(wage_gap))
		#info("Min: ", round(sum(epsilon_per_L[:] .== parameters[:numerical_zero])/length(epsilon_per_L[:])*100,2), "%, ", "Max: ", round(sum(epsilon_per_L[:] .== 1 - parameters[:numerical_zero])/length(epsilon_per_L[:])*100,2), "%")
		
		L_njs = (1 - step_size) * L_njs .+ step_size * (L_nj_star .+ L_n .* epsilon_per_L)
		# make sure no negative L_njs, but they sum to L_n
		L_njs = L_n .* L_njs ./ sum(L_njs, 3)
		#if sum(L_njs[:] .<= 0) > 0
		#	info(L_njs)
		#end

		dist = distance(L_njs, random_variables[:L_njs])
		info("---- Adjustment ", k, ": ", dist)

		k = k+1
	end
	#warn("adjustmenet: ", k-1)
	debug("-- END Adjustment loop")
end

function expected_wage_share(random_variables)
	w_njs = random_variables[:w_njs]
	L_njs = random_variables[:L_njs]

	wage_bill = w_njs .* L_njs
	wage_share = wage_bill ./ sum(wage_bill, 3)

	return expected_value(wage_share)
end

function outer_loop!(random_variables, parameters, t)
	N, J = parameters[:N], parameters[:J]
	lambda = parameters[:outer_step_size]
	L_nj_star = ones(1,N,J,1) / J
	random_variables[:L_njs] = L_nj_star
	starting_values!(random_variables, parameters, t)

	debug("BEGIN Outer loop")
	dist = 999
	k = 1

	while (dist > parameters[:outer_tolerance]) && (k <= parameters[:max_iter_outer])
		old_wage_share = L_nj_star ./ sum(L_nj_star, 3)
		adjustment_loop!(random_variables, L_nj_star, parameters, t)
		wage_share = expected_wage_share(random_variables)

		dist = distance(wage_share, old_wage_share)
		info("-- Outer ", k, ": ", dist)

		L_nj_star = ((1-lambda)*old_wage_share + lambda*wage_share)
		k = k+1
		debug(sum(random_variables[:R_njs][:,:,:,1]))
	end
	#warn("outer: ", k-1)
	debug("END Outer loop")
	return L_nj_star
end

function period_wrapper(parameters, t)
	info("--- Period ", t, " ---")
	A_njs = parameters[:A_njs][t]
	random_variables = Dict{Symbol, Any}()
	random_variables[:A_njs] = A_njs 
	L_nj_star = outer_loop!(random_variables, parameters, t)
	compute_real_gdp!(random_variables, parameters, t)
	return random_variables
end

end