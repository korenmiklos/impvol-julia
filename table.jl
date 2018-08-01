include("output.jl")
parameters = Dict{Symbol, Any}()
parameters[:N] = 25
parameters[:bp_weights] = [0.774074394803123, -0.201004684236153, -0.135080548288772,-0.050951964876636]
stats = ImpvolOutput.write_results(parameters, :real_GDP, true, ARGS[1])
