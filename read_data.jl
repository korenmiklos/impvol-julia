function read_data(relativepath::String, dim_size::NTuple, in_pos::Array{Bool,1}, out_pos::Array{Int64,1}, delim::Char, header::Bool, drop::Int64, indic::Bool, convfloat::Bool)
	# relativepath - gives the relative path of the data to be imported
	# dim_size     - gives the size of each dimension the data should be stored in julia,
	#                their product should correspond to the number of datapoints in the dataset
	#                to be imported, order should be appropriate
	# in_pos       - gives existing dimensions in the input data matrix 
	#                corresponding to mnjt ordering notation as expected (but not restricted) in the output,
	#                its length nevertheless defines the dimensionality of the output matrix
	# out_pos      - gives the position of each dimension in the output matrix
	#                as in the order of argument 'dim_size',
	#                the order prefarabely corresponds to mnjt
	#                m - importer country, n - exporter country, j - sector, t - period
	# delim        - delimiter of the CSV.read command
	# header       - header of the CSV.read command
	# drop         - drop the first 'drop' number of columns in the input data matrix
	# indic        - indicator for reading import shares
	#
	# example: read_data("data/raw_imputed/beta_panel.txt",(36,25,24),[false,true,true,true],[2,3,1],'\t',true,2,false)

	length(dim_size) == count(in_pos) || error("Each dimension should be present in the input matrix")
	any(x -> x <= length(in_pos), out_pos) || error("The output matrix cannot be more than $(length(in_pos))-dimensional")
	unique(out_pos) == out_pos || error("Each output dimension should be unique")

	length(dim_size) == length(out_pos) || error("Each dimension should have its own output position and vice versa")

	data = CSV.read(relativepath, header = header, delim = delim)

	# Part purely for manipulating import shares
	if indic
		data, dim_size = manipulate_import_shares(data,dim_size)
	end
	# End of manipulation

	data = convert(Array,data)
	if drop > 0
		data = data[:,(drop + 1):end]
	end
	prod(size(data)) == prod(dim_size) || error("The number of datapoints ($(prod(size(data)))) should match the product of elements in argument 'dim_size' ($(prod(dim_size)))")

	data = reshape(data, dim_size)

	pos = zeros(Int64,size(in_pos))
	pos[in_pos] = out_pos
	if length(in_pos) > length(dim_size)
		pos[.!in_pos] = convert(Array{Int64,1},length(in_pos):-1:(length(dim_size) + 1))
	end
	pos = convert(Array{Int64,1},pos)
	data = permutedims(cat(ndims(data) + (length(in_pos) - length(dim_size)), data), pos)

	if convfloat
		data = collect(Missings.replace(data,NaN))
		data = convert(Array{Float64,length(in_pos)},data)
	end

	return data
end