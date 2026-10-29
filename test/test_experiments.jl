# Sweep-harness tests: every alternative predicts on the scale it is scored on, and
# coverage is measured against the interval that should contain an observed well.

using Test, LabQC, Statistics, Random
using LabQC: absolute, relative, truth_surface, error_target, baseline_methods,
             Surface, sweep_cell, layout_geometry

const RE, CE = 8, 12
surface_data() = [0.3 + 0.4 * exp(-((r - 4.5)^2 / 10 + (c - 9.0)^2 / 20)) + 0.01 * sin(3r + c)
                  for r in 1:RE, c in 1:CE]

function some_mask(n; seed = 3)
    Random.seed!(seed)
    m = falses(RE, CE)
    m[shuffle(vec(CartesianIndices(m)))[1:n]] .= true
    return m
end

@testset "alternatives predict on the surface's own scale" begin
    G = surface_data(); mask = some_mask(12)
    # the target has no default scale, so the old absolute-only behaviour cannot return
    @test_throws MethodError error_target(G, mask)
    for corr in (absolute, relative)
        t = truth_surface(G, mask, corr)
        for (name, f) in baseline_methods
            P = f(G, mask, corr)
            @test size(P) == size(G)
            # nearest-control and inverse-distance reproduce the target exactly at controls
            name in ("nearest", "idw") && @test P[mask] ≈ t[mask]
        end
    end
    # relative errors are fractions of the well, not OD differences
    t_rel = truth_surface(G, mask, relative)
    P_rel = baseline_methods["idw"](G, mask, relative)
    @test P_rel[mask] ≈ (G[mask] .- median(G[mask])) ./ G[mask]
    @test !(P_rel[mask] ≈ (G[mask] .- median(G[mask])))
end

@testset "coverage uses the predictive interval" begin
    G = surface_data(); mask = some_mask(16)
    s = Surface("synthetic", G, absolute, :negatives, Dict{String,Any}())
    geom = layout_geometry(mask, RE, CE)
    scores = (minisum = NaN, lhs = NaN, hybrid = NaN)
    rows = sweep_cell(s, mask, geom, scores)
    gp = only(filter(r -> r.model == "gp", rows))
    @test gp.mean_sd >= gp.mean_sd_f            # predictive sd includes the noise term
    @test gp.cover95 >= gp.cover95_f
    @test gp.cover68 <= gp.cover95
    @test gp.n_heldout == count(.!mask)
    # with no near-zero wells, the screened set is every held-out well
    @test gp.n_screened == gp.n_heldout
    @test gp.skill_screened ≈ gp.skill
    @test length(unique(map(keys, rows))) == 1   # one schema for every model's row
end
