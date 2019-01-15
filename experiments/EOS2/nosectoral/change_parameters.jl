# parameters that govern counterfactual
parameters[:sigma] = 2.0

remove_shock!(parameters, :global_sectoral_shock_njs)
remove_shock!(parameters, :idiosyncratic_shock_njs)
