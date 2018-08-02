# parameters that govern counterfactual

remove_shock!(parameters, :global_sectoral_shock_njs)
remove_shock!(parameters, :idiosyncratic_shock_njs)

## No labor adjustment
parameters[:one_over_rho] = 0.0