using Test

@testset "LabQC" begin
    include("test_ingest.jl")
    include("test_gp.jl")
end
