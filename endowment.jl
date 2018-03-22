function step(L_nt,alpha_nt)
	return L_nt .^ (-1/theta) .* permutedims(sum(alpha_nt .* tau_nmt .^ (1-theta), 1),[2 1 3]) .^ (1/theta)
end

function Phi(w_nt)
	return sum(tau_nmt .* permutedims(w_nt, [2 1 3]), 2) .^ (1-theta)
end

function alpha(w_nt, L_nt)
	a = w_nt .* L_nt ./ Phi(w_nt)
	return a ./ sum(a, 1)
end

function loop(L_nt)
	alpha_guess = ones(L_nt)/N
	for k in 1:6
		new_wage = step(L_nt, alpha_guess)
		alpha_guess = alpha(new_wage, L_nt)
		println(new_wage[1], "---", new_wage[end])
	end
end

global N=60
global T=10
global theta = 4

global tau_nmt = exp(2.0*rand(N,N,T))
for n=1:N
	for t=1:T
		tau_nmt[n,n,t] = 1.0
	end
end

L_nt = ones(N,N,T) .* [i^3 for i in 1:N]'

@time loop(L_nt)

# free trade
tau_nmt = ones(N,N,T)
@time loop(L_nt)

# autarky
tau_nmt = 100*ones(N,N,T)
for n=1:N
	for t=1:T
		tau_nmt[n,n,t] = 1.0
	end
end
@time loop(L_nt)
