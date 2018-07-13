using CSV
using DataFrames
include("read_data.jl")
include("calibration_utils.jl")
include("calibration_params.jl")

country_names   = read_data("data/raw_imputed/country_name.txt",(25,),[false,true,false,false],[1],'\t',false,0,false,false)
beta            = read_data("data/raw_imputed/beta_panel.txt",(36,25,24),[false,true,true,true],[2,3,1],'\t',true,2,false,true)                
pwt             = read_data("data/raw_imputed/aggregate_price_relative_to_US.csv",(36,25),[false,true,false,true],[2,1],',',false,2,false,true)
va              = read_data("data/raw_imputed/sectoral_value_added.csv",(36,25,24),[false,true,true,true],[2,3,1],',',true,2,false,true)
import_shares   = read_data("data/raw_imputed/import_share.txt",(24,25,36,23),[true,true,true,true],[2,1,4,3],'\t',false,3,true,true)
io_values       = read_data("data/raw_imputed/oecd_io_values.csv",(34,34,13),[true,true,false,true],[2,1,3],',',true,3,false,true)
total_output    = read_data("data/raw_imputed/oecd_total_output.csv",(34,13),[false,true,false,true],[1,2],',',true,2,false,true)
output_shares   = read_data("data/raw_imputed/output_shares.csv",(13,10),[false,true,false,true],[2,1],',',true,1,false,true)
intermediate_input_shares = read_data("data/raw_imputed/intermediate_input_shares.csv",(13,10),[false,true,false,true],[2,1],',',true,1,false,true)
trade_balance   = read_data("data/raw_imputed/trade_balance_new.csv",(25,36),[false,true,false,true],[1,2],',',true,1,false,true)
p_sectoral_data = read_data("data/raw_imputed/sectoral_price_index.csv",(36,18,24),[false,true,true,true],[2,3,1],',',false,0,false,true)

weights = [0.774074394803123; -0.201004684236153; -0.135080548288772; -0.0509519648766360]
io_links = true

if io_links
	gammas = compute_gammas(beta,io_values,total_output,output_shares,intermediate_input_shares)
	alphas = compute_alphas(va,beta,gammas,weights)
else
	beta = ones(size(beta))
	gammas = compute_gammas(beta,io_values,total_output,output_shares,intermediate_input_shares)
	alphas = compute_alphas(va,beta,gammas,weights)
end

d = expenditure_shares(import_shares, 0.000001)

kappa = trade_costs(d, 4.0, 0.000001)



# # Pkg.add("DataFrames")
# # Pkg.add("CSV")

# # using DataFrames
# using CSV

# @everywhere include("utils.jl")

# country_names = CSV.read("data/raw_imputed/country_name.txt", header = false)
# # 25×1, country x 1
# N = size(country_names,1)

# beta = convert(Array, CSV.read("data/raw_imputed/beta_panel.txt", header = true, delim = '\t'))
# beta = beta[:,3:end]
# J = size(beta,2)
# T = Int(size(beta,1)/N)
# beta = reshape(beta,T,N,J) # time,country,sect -> country(25),sect(24),time(36)
# beta = permutedims(cat(ndims(beta) + 1,beta),[4,2,3,1])
# # 25×24×36, country x sect x time
# beta = convert(Array{Float64,4},beta)
# beta = squeeze(mean(beta,(1,2,4)),(1,2,4))
# beta = beta[:,:] # it makes sure, that beta is a matrix and not a vector

# pwt = convert(Array,CSV.read("data/raw_imputed/aggregate_price_relative_to_US.csv", header = false, delim = ','))
# pwt = pwt[:,3]
# pwt = reshape(pwt,T,N) # time,country -> country(25),time(36)
# pwt = permutedims(pwt,[2,1])
# # 25×36, country x time
# pwt = convert(Array{Float64,2},pwt)

# va = convert(Array,CSV.read("data/raw_imputed/sectoral_value_added.csv", header = true, delim = ','))
# va = va[:,3:end]
# va = reshape(va,T,N,J) # time,country,sect -> country(25),sect(24),time(36)
# va = permutedims(va,[2,3,1])
# # 25×24×36, country x sect x time
# va = convert(Array{Float64,3},va)

# import_shares = CSV.read("data/raw_imputed/import_share.txt", header = false, delim = '\t')
# for t in 1:T
# 	y = 1971 + t
# 	for d in 1:N
# 		vec = zeros(1,J - 1 + 3)
# 		vec[1,1:3] = [y,d,d]
# 		push!(import_shares,vec)
# 	end
# end
# import_shares = convert(Array,sort!(import_shares, (1,2,3)))
# import_shares = import_shares[:,4:end]
# import_shares = reshape(import_shares,N,N,T,J - 1) # orig(25),dest(25),time(36),sect(23) -> dest(25),orig(25),sect(23),time(36)
# import_shares = permutedims(import_shares,[2,1,4,3])
# # 25×25×23×36, dest x orig x sect x time
# import_shares = convert(Array{Float64,4},import_shares)

# io_values = convert(Array,CSV.read("data/raw_imputed/oecd_io_values.csv", header = true, delim = ','))
# io_values = io_values[:,4]
# io_values = reshape(io_values,34,34,13) # orig(34),dest(34),time(13) -> dest(34),orig(34),time(13)
# io_values = permutedims(io_values,[2,1,3])
# # 34×34×13, dest_sect x orig_sect x time
# io_values = convert(Array{Float64,3},io_values)

# total_output = convert(Array,CSV.read("data/raw_imputed/oecd_total_output.csv", header = true, delim = ','))
# total_output = total_output[:,3]
# total_output = reshape(total_output,34,13) # dest(34),time(13)
# # 34×13, sector x time
# total_output = convert(Array{Float64,2},total_output)

# output_shares = convert(Array,CSV.read("data/raw_imputed/output_shares.csv", header = true, delim = ','))
# output_shares = output_shares[:,2:end]
# output_shares = permutedims(output_shares,[2,1])
# # 10×13, sector x time
# output_shares = convert(Array{Float64,2},output_shares)

# intermediate_input_shares = convert(Array,CSV.read("data/raw_imputed/intermediate_input_shares.csv", header = true, delim = ','))
# intermediate_input_shares = intermediate_input_shares[:,2:end]
# intermediate_input_shares = permutedims(intermediate_input_shares,[2,1])
# # 10×13, sector x time
# intermediate_input_shares = convert(Array{Float64,2},intermediate_input_shares)

# trade_balance = convert(Array,CSV.read("data/raw_imputed/trade_balance_new.csv", header = true, delim = ','))
# trade_balance = trade_balance[:,2:end]
# # 25×36, country x time
# trade_balance = convert(Array{Float64,2},trade_balance)

# p_sectoral_data = convert(Array,CSV.read("data/raw_imputed/sectoral_price_index.csv", header = false, delim = ','))
# p_sectoral_data = reshape(p_sectoral_data,36,18,24) # time(36),country(18),sectors(24) -> country(18),sectors(24),time(36)
# p_sectoral_data = permutedims(p_sectoral_data,[2,3,1])
# # 18×24×36, country x sector x time

# # IO links
# #if io_links = 1
# #	gammas = compute_gammas(beta,io_values,total_output,output_shares,intermediate_input_shares)
# 	#alpha = 
# 	#gammas = repmat...
# #else

# #end