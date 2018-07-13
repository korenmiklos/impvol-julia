module Scenario
	import ..Environment: N, J, T

	export parameters
	parameters = Dict{Symbol, Any}()
	# parameters that govern counterfactual
	parameters[:kappa_mnjt] = rand(N,N,J,T)
	for j=1:J
		for t=1:T
			for n=1:N
				parameters[:kappa_mnjt][n,n,j,t] = 1.0
			end
		end
	end
end
