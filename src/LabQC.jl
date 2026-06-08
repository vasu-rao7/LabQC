module LabQC

using PlateArrays, DataFrames, CategoricalArrays, GaussianProcesses, Plots, Statistics, StatsBase, Random

include("utils.jl")
include("correction/plateerrorGP.jl")
include("correction/correction.jl")

include("plotting/plotting.jl")
include("reports/control_report.jl")


export data_correction
export plot_control_data, plot_plate_data 
export control_summary

end # module LabQC
