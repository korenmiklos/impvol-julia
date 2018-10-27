.PHONY: data install tables
.PRECIOUS: experiments/*/*/results.jld2
UTILS = calibrate_params.jl calibration_utils.jl utils.jl equilibrium.jl experiments/config.jl
COLUMNS = actual kappa1972 nosectoral nosectoral_kappa1972
TABLES = baseline CES china_1972 no_china no_io_linkages no_labor_adjustment trade_imbalance
PROCS = -p10

tables: $(foreach table,$(TABLES),experiments/$(table)/output_table.csv) 

define run_experiment
experiments/$(1)/%/results.jld2: experiments/$(1)/init_parameters.jl experiments/$(1)/%/scenario.jl experiments/$(1)/%/change_parameters.jl $(UTILS)  data/impvol_data.jld2 
	@echo " + Compiling '$$@'"
	cd $$(dir $$@) && julia $(PROCS) scenario.jl
endef

$(foreach experiment,$(TABLES),$(eval $(call run_experiment,$(experiment))))

experiments/%/output_table.csv: $(foreach column,$(COLUMNS),experiments/%/$(column)/results.jld2) output.jl table.jl
	julia table.jl $(dir $@)

data: data/impvol_data.jld2
data/impvol_data.jld2: read_data.jl data/*.csv data/*.txt
	julia read_data.jl

install: install.jl
	julia install.jl
