.PHONY: data install test
UTILS = calibrate_params.jl calibration_utils.jl utils.jl equilibrium.jl

test: experiments/test/scenario1/results.jld

experiments/test/scenario1/results.jld: $(UTILS)  data/impvol_data.jld experiments/test/config.jl experiments/test/scenario1/parameters.jl experiments/test/scenario1/scenario.jl
	cd experiments/test/scenario1; julia scenario.jl

data: data/impvol_data.jld
data/impvol_data.jld: read_data.jl data/*.csv data/*.txt
	julia read_data.jl

install: install.jl
	julia install.jl