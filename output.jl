module ImpvolOutput
	using Plots
	include("calibration_utils.jl")

	function calculate_volatilities(x::Array{Float64,2}, parameters, bool_detrend::Bool)
		if bool_detrend
			x_c, x_t = detrend(log(x),parameters[:bp_weights])
		else
			x_c = log(x)
		end

		return var(x_c,2)
	end

	function plot_model_vs_data(data::Array{Float64,2}, model::Array{Float64,2}, label, title::String)
		# data and model are expected to be of dimension N x T
		# length of label should be equal N
		
		data = transpose(data)
		model = transpose(model)

		size(data,2) == size(model,2) || error("You can only compare matching time series between the model outcome and data")

		colors = distinguishable_colors(size(data,2))

		fig = plot()
		for i in 1:length(colors)
			plot!([data[:,i] model[:,i]], color = colors[i], ls = [:solid :dash], label = [label[i] ""], title = title)
		end
		return fig
	end
end