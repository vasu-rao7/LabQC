# The placement / control-count sweep.
#
# One loop produces all three manuscript experiments, because they read the same
# fits: layout comparison is the sweep at a fixed control count, the accuracy curve
# is the sweep across control counts, and the prediction-vs-truth scores come from
# every fit. Results stream to CSV so a long run survives interruption, and a rerun
# skips rows that already exist.

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

"""
    sweep_cell(surface, mask, geom, scores; pitch, ls, lb_ls, ub_ls) -> Vector{NamedTuple}

Fit the GP and every baseline to one layout on one surface, scoring all of them on
the held-out wells. Control wells are excluded from scoring: the GP is pinned there
by construction, so including them inflates every model identically and hides the
differences the sweep exists to measure.
"""
function sweep_cell(s::Surface, mask::BitMatrix, geom, scores;
                    pitch = well_pitch(size(s.values)...),
                    ls = 3 * pitch, lb_ls = 2 * pitch, ub_ls = 6 * pitch,
                    screen_frac = 0.25)
    truth = truth_surface(s.values, mask, s.correction)
    held = .!mask
    t = truth[held]
    # wells that barely grew make relative error explode (0.04 vs 0.55 is about -12),
    # so every score is also reported on held-out wells at or above a fraction of the
    # surface median; on surfaces without such wells this is the full held-out set
    screened = held .& (s.values .>= screen_frac * median(s.values))
    ts = truth[screened]
    rows = NamedTuple[]

    el = @elapsed begin
        _, diag = plateerrorGP(s.values, mask, s.correction;
                               verbose = true, pitch = pitch, ls = ls,
                               lb_ls = lb_ls, ub_ls = ub_ls)
    end
    gp = diag["gp"]
    preds = Dict{String,Matrix{Float64}}("gp" => diag["mean_error"])
    for (name, f) in baseline_methods
        preds[name] = f(s.values, mask, s.correction)
    end

    # Held-out truth is an observed deviation, so it carries the per-well noise the GP
    # fitted. The interval that should contain it is the predictive one (var_y); the
    # latent one (var_f) covers only the smooth surface and is kept for comparison.
    sd_y = sqrt.(diag["var_y"])
    sd_f = diag["sd_error"]
    ell = lengthscales(gp)
    nan = fill(NaN, count(held))
    for (model, P) in preds
        isgp = model == "gp"
        p, ps = P[held], P[screened]
        sy, sf = isgp ? sd_y[held] : nan, isgp ? sd_f[held] : nan
        push!(rows, (model = model,
                     skill = skill(p, t), rmse = rmse(p, t), mae = mae(p, t),
                     bias = bias(p, t), spearman = spearman(p, t),
                     cover68 = isgp ? coverage(p, t, sy; z = 1.0) : NaN,
                     cover95 = isgp ? coverage(p, t, sy; z = 1.96) : NaN,
                     mean_sd = isgp ? mean(sy) : NaN,
                     max_sd = isgp ? maximum(sy) : NaN,
                     cover68_f = isgp ? coverage(p, t, sf; z = 1.0) : NaN,
                     cover95_f = isgp ? coverage(p, t, sf; z = 1.96) : NaN,
                     mean_sd_f = isgp ? mean(sf) : NaN,
                     max_sd_f = isgp ? maximum(sf) : NaN,
                     skill_screened = skill(ps, ts),
                     cover95_screened = isgp ? coverage(ps, ts, sd_y[screened]; z = 1.96) : NaN,
                     n_screened = count(screened),
                     ls_fit = ell[1], sd_fit = signal_sd(gp),
                     noise_fit = exp(gp.logNoise.value),
                     max_min_dist = geom.max_min_dist, mean_min_dist = geom.mean_min_dist,
                     cols_covered = geom.cols_covered, rows_covered = geom.rows_covered,
                     minisum = scores.minisum, lhs_score = scores.lhs,
                     hybrid_score = scores.hybrid,
                     n_heldout = count(held), seconds = el))
    end
    return rows
end
