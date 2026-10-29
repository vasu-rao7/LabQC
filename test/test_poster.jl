# Poster plate E4: the 2021 correction must reproduce, and the replicate score must be
# blind to shifts and rescaling so that only position-dependent changes move it.

using Test, LabQC, Statistics
using LabQC: read_mantis_layout, two_pass, replicate_scores, poster_metrics, control_subsets

const E4 = joinpath(@__DIR__, "..", "raw_data", "poster_e4")

@testset "poster plate E4" begin
    k = read_kinetic(joinpath(E4, "E4_211118_CDM_Glucose_V2_cleaned.xlsx"))
    lay = read_mantis_layout(joinpath(E4, "E4_layout_SMU_Glucose_Mantis_Input2_111621.csv"))
    @test size(lay.conc) == (16, 24)
    @test count(lay.negatives) == 24 && count(lay.positives) == 24 && count(lay.middle) == 336
    @test length(unique(lay.conc[lay.middle])) == 49
    m, _ = poster_metrics(k)
    G = m["od_read90"]

    # the 2021 notebook's printed outputs (cells 56 and 62), ARD kernel, read 90
    X, info = two_pass(G, lay.negatives, lay.positives, "gp_ard")
    @test info["p1_mean_error"][1, 1] ≈ -0.00017611731471516012 atol = 1e-8
    @test info["p1_mean_error"][16, 24] ≈ -0.010694503026647761 atol = 1e-8
    @test X[1, 1] ≈ 0.15733070657529374 atol = 1e-8
    @test X[16, 24] ≈ 0.32016871294467675 atol = 1e-8
    @test info["p2_logparams"] ≈ [1.7917594692213736, 0.693147180614958, -2.927836001984723] atol = 1e-6

    # replicate skill ignores shifts and rescaling: raw and any affine normalisation score 0
    s_raw = replicate_scores(G, G, lay.conc, lay.middle, lay.negatives, lay.positives)
    @test s_raw.replicate_skill ≈ 0 atol = 1e-12
    Y, _ = two_pass(G, lay.negatives, lay.positives, "subtract_and_normalize")
    s_sn = replicate_scores(Y, G, lay.conc, lay.middle, lay.negatives, lay.positives)
    @test s_sn.replicate_skill ≈ 0 atol = 1e-12
    @test s_sn.eta2 ≈ s_raw.eta2
    @test s_sn.within_sd ≈ s_raw.within_sd     # window units are scale-free too
    # perfect agreement scores 1
    perfect = [lay.middle[i, j] ? mean(G[lay.middle .& (lay.conc .== lay.conc[i, j])]) : G[i, j]
               for i in 1:16, j in 1:24]
    @test replicate_scores(perfect, G, lay.conc, lay.middle, lay.negatives, lay.positives).replicate_skill ≈ 1

    # subsets: paired, the right sizes, drawn from the plate's own controls
    subs = control_subsets(lay.negatives, lay.positives, [3, 24], 5)
    @test length(subs) == 5 + 1                  # k = 24 is the full set, drawn once
    @test all(s -> count(s.negatives) == s.k && count(s.positives) == s.k, subs)
    @test all(s -> all(lay.negatives[s.negatives]) && all(lay.positives[s.positives]), subs)
    @test subs == control_subsets(lay.negatives, lay.positives, [3, 24], 5)   # reproducible
end
