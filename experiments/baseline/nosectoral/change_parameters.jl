# parameters that govern counterfactual

Logging.debug(size(parameters[:global_sectoral_shock]))

remove_shock!(parameters, :global_sectoral_shock_njs)
remove_shock!(parameters, :idiosyncratic_shock_njs)
