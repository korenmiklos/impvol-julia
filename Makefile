.PHONY: data install tables template calibrate
CALIBRATION = calibrate_params.jl calibration_utils.jl experiments/config.jl data/impvol_data.jld2
EQULIBRIUM = utils.jl equilibrium.jl experiments/config.jl
COLUMNS = actual kappa1972 nosectoral nosectoral_kappa1972
TABLES = baseline CES05 CES2 china_1972 no_china no_io_linkages labor_adjustment trade_imbalance S1000 theta2 theta8
.PRECIOUS: $(foreach table,$(TABLES),$(foreach column,$(COLUMNS),experiments/$(table)/$(column)/results.jld2))
PROCS = -p10

tables: $(foreach table,$(TABLES),experiments/$(table)/output_table.csv) 

calibrate: $(foreach table,$(TABLES),experiments/$(table)/common_parameters.jld2) 

experiments/%/common_parameters.jld2: experiments/%/init_parameters.jl $(CALIBRATION) 
	cd $(dir $@) && julia init_parameters.jl

define run_experiment
experiments/$(1)/%/results.jld2: $(EQULIBRIUM) experiments/$(1)/common_parameters.jld2 experiments/$(1)/%/scenario.jl experiments/$(1)/%/change_parameters.jl 
	@echo " + Compiling '$$@'"
	cd $$(dir $$@) && julia $(PROCS) scenario.jl
endef

$(foreach experiment,$(TABLES),$(eval $(call run_experiment,$(experiment))))

experiments/%/output_table.csv: $(foreach column,$(COLUMNS),experiments/%/$(column)/results.jld2) output.jl table.jl
	julia table.jl $(dir $@)

data: data/impvol_data.jld2
data/impvol_data.jld2: read_data.jl data/*.csv data/*.txt
	julia read_data.jl

template: scenario_template.jl
	find . -name "scenario.jl" -exec cp scenario_template.jl {} \; 

install: install.jl
	julia install.jl
