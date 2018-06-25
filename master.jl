@everywhere using Logging

@everywhere include("utils.jl")
@everywhere include("equilibrium.jl")
@everywhere using ImpvolEquilibrium

parameters = Dict{Symbol, Any}()

# for matrix conformity, store all variables in a 4-dimensional array:
# mnjt: destination, source, sector, time

# per-period random variables are stored as
# mnjs: destination, source, sector, state


@everywhere Logging.configure(level=DEBUG)

N = 2
J = 3
T = 10
S = 1000
# set random seed for reproducability
srand(2311)
fill_dict!(parameters, N=N, J=J, T=T, S=S)

alpha = rand(J) .* ones(J,T)
beta = 0.25 + 0.75 * rand(1, J)
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
# make Gamma more diagonal
gamma_jk = rand(J,J) + 1.0*eye(J)
# for testing purposed, set IO links to 0
# test for continuity with small IO links
#parameters[:gamma_jk] = 0.75*eye(J)
#parameters[:gamma_jk] = repmat((1-beta[:]')/J, J, 1)
#parameters[:beta] = 0.25*ones(1,J)
parameters[:gamma_jk] = gamma_jk ./ sum(gamma_jk, 1) .* (1-beta)
# adaptive step size. large lambda means large steps
parameters[:lambda] = exp(-0.05*(J-1)^0.75)
# inverse of adjustment cost, 0 if cannot readjust
parameters[:one_over_rho] = 0.01
# this is log points of average input price differences
parameters[:inner_tolerance] = 0.001
parameters[:middle_tolerance] = 0.003
parameters[:adjustment_tolerance] = 0.003
parameters[:outer_tolerance] = 0.005
parameters[:numerical_zero] = 1e-6

# standard deviation for each (n,j)
parameters[:sigma] = 0.1*ones(1,N,J,1)
# AR coefficient for each (n,j)
parameters[:AR_decay] = 0.9*ones(1,N,J,1)

coerce_parameters!(parameters)

A_njs = 1.0 .+ rand(1,N,J,S)

@time results = pmap(t -> period_wrapper(A_njs, parameters, t), 1:T)

#for t = 1:T
#	@time random_variables = period_wrapper(A_njs, parameters, t)
#	real_GDP = non_random_variable(random_variables[:real_GDP], 1)
#	A_njs = draw_next_productivity(A_njs, parameters, 1)
#end