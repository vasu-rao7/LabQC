# Control-layout generation for the placement sweep.
#
# Upstream caveats this module works around (all verified against PlateArrays
# v0.1.1, and none of them patched here):
#
#   * `exchange` picks its best restart with `findmax` regardless of the `minimize`
#     flag, so with the default `restarts=10` it returns the WORST of ten local
#     optima and gets worse as `restarts` grows. Every call here uses `restarts=1`,
#     where the bug cannot bite, and draws many independent layouts instead.
#   * `hybrid`'s normalisation bounds are recomputed by four randomised solves on
#     every call and never returned, so lambda is not comparable between calls.
#     `fixed_hybrid_bounds` computes them once so scores share a scale.
#   * `objectives["minimax"]` is a Manhattan-distance *minisum*, not a minimax.
#     `layout_geometry` records true Euclidean descriptors alongside it.

"""
    fixed_hybrid_bounds(wells, P, N; reps=8, iterations=1000, seed=0) -> NamedTuple

Normalisation bounds for the hybrid objective, computed once and reused for the
whole sweep so that a given lambda means the same thing across every call.
"""
function fixed_hybrid_bounds(wells::BitMatrix, P::Int, N::Int;
                             reps = 8, iterations = 1000, seed = 0)
    Random.seed!(seed)
    lo_d, hi_d, lo_l, hi_l = Inf, -Inf, Inf, -Inf
    for _ in 1:reps
        for (obj, minim) in ((PlateArrays.minimax, true), (PlateArrays.minimax, false),
                             (PlateArrays.LHS, true), (PlateArrays.LHS, false))
            pa = PlateArrays.exchange(wells, P, N; objective = obj, minimize = minim,
                                      restarts = 1, iterations = iterations)
            d, l = PlateArrays.minimax(pa), PlateArrays.LHS(pa)
            lo_d = min(lo_d, d); hi_d = max(hi_d, d)
            lo_l = min(lo_l, l); hi_l = max(hi_l, l)
        end
    end
    return (lb_dist = lo_d, ub_dist = hi_d, lb_LHS = lo_l, ub_LHS = hi_l)
end

"""
    generate_layouts(wells, P, N; objective, lambda, nlayouts, seed, iterations)

`nlayouts` independent control layouts from one objective. Each is a separate
single-restart solve, so the collection is a sample of what that objective produces
rather than a single point estimate — which is what makes the layout, not the
layout-timepoint pair, the unit of replication downstream.

`objective` is `"minimax"`, `"LHS"`, `"hybrid"`, or `"random"` for the null.
"""
function generate_layouts(wells::BitMatrix, P::Int, N::Int;
                          objective::String = "hybrid", lambda = 0.5,
                          nlayouts::Int = 50, seed::Int = 1, iterations::Int = 1000)
    out = PlateArray[]
    for k in 1:nlayouts
        Random.seed!(hash((seed, objective, lambda, k)) % typemax(Int32))
        pa = if objective == "random"
            PlateArrays.random_platearray(wells, P, N)
        elseif objective == "hybrid"
            place_controls(wells, P, N; objective = "hybrid", lambda = lambda,
                           restarts = 1, iterations = iterations)
        else
            place_controls(wells, P, N; objective = objective,
                           restarts = 1, iterations = iterations)
        end
        push!(out, pa)
    end
    return out
end

"""
    layout_geometry(mask, nrow, ncol; pitch=well_pitch(nrow,ncol)) -> NamedTuple

Descriptors of one control class's positions, computed here rather than taken from
PlateArrays so that results stay interpretable regardless of upstream naming. All
distances are true Euclidean, in mm.
"""
function layout_geometry(mask::BitMatrix, nrow::Int, ncol::Int;
                         pitch = well_pitch(nrow, ncol))
    idx = findall(mask)
    isempty(idx) && return (max_min_dist = NaN, mean_min_dist = NaN,
                            rows_covered = 0, cols_covered = 0, max_per_row = 0, max_per_col = 0)
    dmins = Float64[]
    for i in 1:nrow, j in 1:ncol
        push!(dmins, minimum(sqrt((ci[1] - i)^2 + (ci[2] - j)^2) * pitch for ci in idx))
    end
    rows = [count(ci -> ci[1] == i, idx) for i in 1:nrow]
    cols = [count(ci -> ci[2] == j, idx) for j in 1:ncol]
    return (max_min_dist = maximum(dmins),      # true minimax (worst-served well)
            mean_min_dist = mean(dmins),        # I-optimality flavour
            rows_covered = count(>(0), rows),
            cols_covered = count(>(0), cols),
            max_per_row = maximum(rows),
            max_per_col = maximum(cols))
end

"""Objective values recomputed on fixed bounds, so they are comparable across calls."""
function layout_scores(pa::PlateArray, bounds; lambda = 0.5)
    return (minisum = PlateArrays.minimax(pa),
            lhs = PlateArrays.LHS(pa),
            hybrid = PlateArrays.hybrid(pa; lambda = lambda, bounds...))
end
