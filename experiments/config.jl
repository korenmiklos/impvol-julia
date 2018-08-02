using Logging
Logging.configure(level=INFO)
include("../calibrate_params.jl")
using .CalibrateParameters

parameters[:S] = 100

parameters[:numerical_zero] = 1e-12

parameters[:bp_weights] = [0.774074394803123; -0.201004684236153; -0.135080548288772; -0.0509519648766360]

CalibrateParameters.calibrate_parameters!(parameters)

# adaptive step size. large lambda means large steps
parameters[:inner_step_size] = exp(-0.10*(parameters[:J]-1)^0.75)
# large substitution needs more dampening
parameters[:middle_step_size] = exp(-0.275*max(1,parameters[:sigma]))
# any deviation from sigma=1 needs more dampening
parameters[:outer_step_size] = exp(-0.5*abs(log(parameters[:sigma])))
# this is log points of average input price differences
parameters[:inner_tolerance] = 0.001
parameters[:middle_tolerance] = 0.003
parameters[:adjustment_tolerance] = 0.003
parameters[:outer_tolerance] = 0.005

