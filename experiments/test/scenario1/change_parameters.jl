# parameters that govern counterfactual
parameters[:kappa_mnjt] = rand(N,N,J,T)
for j=1:J
	for t=1:T
		for n=1:N
			parameters[:kappa_mnjt][n,n,j,t] = 1.0
		end
	end
end
A_njs = 1.0 .+ rand(1,N,J,S)
# to test CES
A_njs[1,1,1,1] = 1.0
A_njs[1,1,2,1] = 0.5
