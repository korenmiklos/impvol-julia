.PHONY: data install table1
.PRECIOUS: experiments/baseline/%/results.jld2
UTILS = calibrate_params.jl calibration_utils.jl utils.jl equilibrium.jl experiments/config.jl
COLUMNS = actual kappa1972 nosectoral nosectoral_kappa1972
PROCS = -p10

table1: experiments/baseline/output_table.csv
experiments/%/output_table.csv: $(foreach column,$(COLUMNS),experiments/%/$(column)/results.jld2) output.jl table.jl
	julia table.jl $(dir $@)

experiments/baseline/%/results.jld2: $(UTILS)  data/impvol_data.jld2 experiments/baseline/init_parameters.jl experiments/baseline/%/scenario.jl experiments/baseline/%/change_parameters.jl
	cd $(dir $@) && julia $(PROCS) scenario.jl

experiments/no_labor_adjustment/%/results.jld2: $(UTILS)  data/impvol_data.jld2 experiments/no_labor_adjustment/init_parameters.jl experiments/no_labor_adjustment/%/scenario.jl experiments/no_labor_adjustment/%/change_parameters.jl
	cd $(dir $@) && julia $(PROCS) scenario.jl

data: data/impvol_data.jld2
data/impvol_data.jld2: read_data.jl data/*.csv data/*.txt
	julia read_data.jl

install: install.jl
	julia install.jl
