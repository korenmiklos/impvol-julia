module ImpvolOutput
	using Plots
	using JLD2
	using FileIO
	using CSV
	using Missings
	include("calibration_utils.jl")

	function read_results(path = "experiments/baseline/actual/results.jld2")
		return load(path)["results"]
	end

	function sort_results(results)
		return sort(collect(results), by = x -> x[1])
	end

	function list_keys(results)
		for (key, value) in results[1][2]
			println(key)
		end
	end

	function make_series(results, key)
		series = zeros(size(results[1][2][key])[1], size(results[1][2][key])[2], size(results[1][2][key])[3], length(results))
		for t in 1:length(results)
			# Only the first element in the shock dimension is interesting, the rest are there only for optimizaton purposes used in the algorithm
			series[:,:,:,t] = results[t][2][key][:,:,:,1]
		end
		return series
	end

	function calculate_volatilities(x::Array{Float64,2}, parameters, bool_detrend::Bool)
		if bool_detrend
			x_c, x_t = detrend(log(x),parameters[:bp_weights])
		else
			x_c = log(x)
		end

		return var(x_c,2)
	end

	function plot_model_vs_data(plot_data::Tuple, title::String)
		# data and model are expected to be of dimension N x T
		# length of label should be equal N
		data = plot_data[1]
		model = plot_data[2]
		label = plot_data[3]

		data = transpose(data)
		model = transpose(model)

		info(size(data,2))
		info(size(model,2))
		size(data,2) == size(model,2) || error("You can only compare matching time series between the model outcome and data")

		colors = distinguishable_colors(size(data,2))

		fig = plot()
		for i in 1:size(data,2)
			plot!([data[:,i] model[:,i]], color = colors[i,1], ls = [:solid :dash], label = [label[i,1] ""], title = title)
		end
		return fig
	end


	function plot_data(key, country_range = ":", path = "experiments/baseline/actual/results.jld2")
		# output: "data, model, label", as an input for the function 'plot_model_vs_data'
		data = load("data/impvol_data.jld2")
		gdp_d = squeeze(sum(data["va"],3), (1,3))
		info(size(gdp_d))
		gdp_d = gdp_d[eval(parse(country_range)),:]
		gdp_d = Float64.(collect(Missings.replace(gdp_d, NaN)))

		cpi = CSV.read("data/cpi.csv", header = true)
		cpi = permutedims(convert(Array, cpi)[:,2:end], (2,1))
		cpi = cpi ./ cat(2, cpi[:,24]) # 24 = 1995-base
		info(size(cpi))
		cpi = cpi[eval(parse(country_range)),:]
		cpi = Float64.(collect(Missings.replace(cpi, NaN)))

		xr = CSV.read("data/exchange_rates.csv", header = true)
		xr = permutedims(convert(Array, xr)[:,2:end], (2,1))
		info(size(xr))
		xr = xr[eval(parse(country_range)),:]
		xr = Float64.(collect(Missings.replace(xr, NaN)))

		country_names = CSV.read("data/country_name.txt", header = false)
		country_names = convert(Array, country_names)
		info(size(country_names))
		country_names = country_names[eval(parse(country_range))]

		gdp_m = make_series(sort_results(read_results(path)), key)
		gdp_m = sum(gdp_m[1,eval(parse(country_range)),:,:],2)
		gdp_m = squeeze(gdp_m, 2)

		return log.(gdp_d ./ cpi .* xr), log.(gdp_m), country_names
	end
end