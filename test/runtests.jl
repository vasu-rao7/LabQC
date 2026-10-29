using Test

@testset "LabQC" begin
    include("test_ingest.jl")
    include("test_gp.jl")
    include("test_experiments.jl")
    include("test_poster.jl")
end
