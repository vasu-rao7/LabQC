# The placement / control-count sweep.
#
# One loop produces all three manuscript experiments, because they read the same
# fits: layout comparison is the sweep at a fixed control count, the accuracy curve
# is the sweep across control counts, and the prediction-vs-truth scores come from
# every fit. Results stream to CSV so a long run survives interruption, and a rerun
# skips rows that already exist.

"""
    truth_surface(G, mask, correction) -> Matrix

The error each model is asked to predict, on the same scale the correction removes:
deviation from the control median for [`absolute`](@ref), and the fractional version
of that for [`relative`](@ref).
"""
function truth_surface(G::AbstractMatrix, mask::BitMatrix, correction::Function)
    A = Float64.(G)
    m = median(A[mask])
    return correction === absolute ? A .- m : (A .- m) ./ A
end

"""
    Surface(id, values, correction, class; meta...)

One measured plate surface to sweep over: a fully-observed all-control plate, the
correction whose scale its errors live on, and which control class a layout should
draw from.
"""
struct Surface
    id::String
    values::Matrix{Float64}
    correction::Function
    class::Symbol            # :negatives or :positives
    meta::Dict{String,Any}
end

const SWEEP_COLUMNS = [:surface, :plate, :metric, :n, :objective, :lambda, :layout,
                       :model, :skill, :rmse, :mae, :bias, :spearman, :cover68, :cover95,
                       :mean_sd, :max_sd, :ls_fit, :sd_fit, :noise_fit,
                       :max_min_dist, :mean_min_dist, :cols_covered, :rows_covered,
                       :minisum, :lhs_score, :hybrid_score, :n_heldout, :seconds]

"""
    sweep_cell(surface, mask, geom, scores; pitch, ls, lb_ls, ub_ls) -> Vector{NamedTuple}

Fit the GP and every baseline to one layout on one surface, scoring all of them on
the held-out wells. Control wells are excluded from scoring: the GP is pinned there
by construction, so including them inflates every model identically and hides the
differences the sweep exists to measure.
"""
function sweep_cell(s::Surface, mask::BitMatrix, geom, scores;
                    pitch = well_pitch(size(s.values)...),
                    ls = 3 * pitch, lb_ls = 2 * pitch, ub_ls = 6 * pitch)
    truth = truth_surface(s.values, mask, s.correction)
    held = .!mask
    t = truth[held]
    rows = NamedTuple[]

    el = @elapsed begin
        _, diag = plateerrorGP(s.values, mask, s.correction;
                               verbose = true, pitch = pitch, ls = ls,
                               lb_ls = lb_ls, ub_ls = ub_ls)
    end
    gp = diag["gp"]
    preds = Dict{String,Matrix{Float64}}("gp" => diag["mean_error"])
    for (name, f) in baseline_methods
        preds[name] = f(s.values, mask)
    end

    for (model, P) in preds
        p = P[held]
        sd = model == "gp" ? diag["sd_error"][held] : fill(NaN, length(p))
        push!(rows, (model = model,
                     skill = skill(p, t), rmse = rmse(p, t), mae = mae(p, t),
                     bias = bias(p, t), spearman = spearman(p, t),
                     cover68 = model == "gp" ? coverage(p, t, sd; z = 1.0) : NaN,
                     cover95 = model == "gp" ? coverage(p, t, sd; z = 1.96) : NaN,
                     mean_sd = model == "gp" ? mean(sd) : NaN,
                     max_sd = model == "gp" ? maximum(sd) : NaN,
                     ls_fit = gp.kernel.ℓ, sd_fit = sqrt(gp.kernel.σ2),
                     noise_fit = exp(gp.logNoise.value),
                     max_min_dist = geom.max_min_dist, mean_min_dist = geom.mean_min_dist,
                     cols_covered = geom.cols_covered, rows_covered = geom.rows_covered,
                     minisum = scores.minisum, lhs_score = scores.lhs,
                     hybrid_score = scores.hybrid,
                     n_heldout = count(held), seconds = el))
    end
    return rows
end
