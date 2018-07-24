.PHONY: data install test
UTILS = calibrate_params.jl calibration_utils.jl utils.jl equilibrium.jl

test: experiments/test/scenario1/results.jld2

experiments/test/scenario1/results.jld2: $(UTILS)  data/impvol_data.jld2 experiments/test/config.jl experiments/test/scenario1/parameters.jl experiments/test/scenario1/scenario.jl
	cd experiments/test/scenario1 && julia scenario.jl

data: data/impvol_data.jld2
data/impvol_data.jld2: read_data.jl data/*.csv data/*.txt
	julia read_data.jl

install: install.jl
	julia install.jl