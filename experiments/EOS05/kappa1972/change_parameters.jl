# parameters that govern counterfactual
parameters[:sigma] = 0.5

## kappa remains at 1972 level
for t=1:parameters[:T]
	parameters[:kappa_mnjt][:,:,:,t] = parameters[:kappa_mnjt][:,:,:,1]
end
