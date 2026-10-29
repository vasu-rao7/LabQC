module LabQC

using PlateArrays, DataFrames, CategoricalArrays, GaussianProcesses, Plots, Statistics, StatsBase, Random
using XLSX, SHA, Dates, LinearAlgebra, CSV, Serialization

include("utils.jl")
include("geometry.jl")
include("io/biotek.jl")
include("correction/plateerrorGP.jl")
include("correction/correction.jl")

include("plotting/plotting.jl")
include("reports/control_report.jl")
include("experiments/metrics.jl")
include("experiments/baselines.jl")
include("experiments/layouts.jl")
include("experiments/sweep.jl")
include("experiments/runner.jl")


export data_correction, gp_correction
export plot_control_data, plot_plate_data 
export control_summary
export read_endpoint_grids, read_kinetic, transient_reads, grid_hash
export well_pitch
export skill, rmse, mae, bias, spearman, coverage, variance_decomposition
export error_target, baseline_methods, lengthscales, signal_sd
export generate_layouts, fixed_hybrid_bounds, layout_geometry, layout_scores
export Surface, truth_surface, sweep_cell
export build_surfaces, default_objectives, layout_bank, run_experiment, EXPERIMENTS

end # module LabQC
