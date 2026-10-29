# Figures for the three experiments, drawn from the tables analyze.jl writes.
#
# Run after experiments/analyze.jl:
#   julia --project experiments/figures.jl
#
# Writes to results/ (PNG for viewing, PDF for the manuscript):
#   fig_exp1_placement      GP skill per layout, by placement objective (n = 12, 40)
#   fig_exp2_coverage       GP 95% interval coverage vs number of controls
#   fig_exp2_sd_vs_error    predictive sd vs realised RMSE, one point per layout
#   fig_exp3_skill_vs_n     skill vs number of controls, GP against simple baselines
#
# Each plate is shown on the scale its correction uses: the negative plate's error is
# additive (technical), the positive and gas plates' error is relative (biological),
# matching `gp_correction`. The gas plates are also shown on the additive scale: on the
# relative scale 2-8 wells per plate that barely grew (relative error near -12) carry
# about 90% of the squared error, both as scored wells and when a control lands on one,
# so the relative panels mostly measure where those wells fall.
#
# Coverage and the sd-vs-error panels use the predictive sd (latent sd plus the
# fitted per-well noise), because held-out truth is an observed well, noise included.

using CSV, DataFrames, Statistics, StatsBase, Plots, Printf, Random
gr()

const OUT = joinpath(@__DIR__, "..", "results")

# (experiment, correction scale, title)
const PANELS = [("negative", "absolute", "All-negative plate (additive)"),
                ("positive", "relative", "All-positive plate (relative)"),
                ("gas",      "relative", "Gas plates (relative)"),
                ("gas",      "absolute", "Gas plates (additive, for comparison)")]
const UNIT = Dict("absolute" => "OD", "relative" => "fraction of well")

# Categorical slots, fixed order (validated: adjacent and 3-way all-pairs CVD-safe; the
# fourth slot is used only for lines, which are also direct-labelled).
const C1, C2, C3, C4 = colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a", colorant"#eda100"
const INK, MUTED, RULE = colorant"#0b0b0b", colorant"#52514e", colorant"#8a8985"

default(fontfamily = "sans-serif", framestyle = :axes, grid = true, gridalpha = 0.12,
        gridlinewidth = 0.5, foreground_color_axis = RULE, foreground_color_border = RULE,
        tickfontcolor = MUTED, guidefontcolor = INK, titlefontcolor = INK,
        titlefontsize = 11, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8,
        legend_foreground_color = nothing, background_color = :white, dpi = 200)

# control-count axes are linear: the sweep runs 3-40, so a log axis stretched the
# few-control end and hid how much is still gained between 24 and 40
const NTICKS = 0:10:40

function save(p, name)
    savefig(p, joinpath(OUT, name * ".png"))
    savefig(p, joinpath(OUT, name * ".pdf"))
    @info "wrote" figure = name
end

# Direct labels at line ends, nudged apart so lines that finish close together stay legible.
function spread(ys, gap)
    order = sortperm(ys)
    out = float.(ys[order])
    for i in 2:length(out)
        out[i] = max(out[i], out[i-1] + gap)
    end
    out .-= (sum(out) - sum(ys)) / length(out)
    res = similar(out); res[order] = out
    return res
end

by_n   = CSV.read(joinpath(OUT, "summary_by_n.csv"), DataFrame)
layout = CSV.read(joinpath(OUT, "per_layout_gp.csv"), DataFrame)
surf   = CSV.read(joinpath(OUT, "surfaces.csv"), DataFrame)
panel(df, e, sc) = filter(r -> r.experiment == e && r.correction == sc, df)

# Ceiling: no smooth interpolator recovers the nugget, so skill is bounded by the
# structured share of variance (median over a plate type's surfaces), on the same scale.
ceiling = Dict((e, sc) => median(panel(surf, e, sc).structured_share) for (e, sc, _) in PANELS)

# --- Experiment 3: skill vs number of controls ---------------------------------
const MODELS = [("gp", "GP (this method)", C1, :circle),
                ("idw", "inverse-distance", C2, :utriangle),
                ("rowcol", "row + column", C3, :rect),
                ("nearest", "nearest control", C4, :diamond)]
const YLIM3 = (-0.7, 0.85)

panels = map(enumerate(PANELS)) do (i, (e, sc, title))
    p = plot(title = title, xticks = NTICKS, xlims = (0, 42), ylims = YLIM3,
             xlabel = "number of control wells", ylabel = i == 1 ? "skill vs. control median" : "",
             legend = i == 1 ? :topleft : false)
    hline!(p, [0], color = RULE, lw = 1, label = "")
    c = ceiling[(e, sc)]
    hline!(p, [c], color = MUTED, lw = 1, ls = :dash, label = "")
    annotate!(p, 0.5, c + 0.035, text(@sprintf("ceiling %.2f", c), 7, MUTED, :left))
    noted = 0
    for (m, lab, col, shape) in MODELS
        d = sort(filter(r -> r.model == m, panel(by_n, e, sc)), :n)
        plot!(p, d.n, d.skill_mean, color = col, lw = 2, marker = shape, ms = 4,
              markerstrokewidth = 0, label = lab)
        if minimum(d.skill_mean) < YLIM3[1]
            annotate!(p, 41.5, YLIM3[2] - 0.05 - 0.07noted,
                      text(@sprintf("%s: min %.1f, off-scale", lab, minimum(d.skill_mean)), 7, MUTED, :right))
            noted += 1
        end
    end
    p
end
save(plot(panels..., layout = (1, 4), size = (1500, 410), left_margin = 4Plots.mm,
          bottom_margin = 5Plots.mm), "fig_exp3_skill_vs_n")

# --- Experiment 1: placement objectives ----------------------------------------
const OBJ = [("random", NaN, "random"), ("LHS", NaN, "LHS"), ("minimax", NaN, "minisum*"),
             ("hybrid", 0.0, "hybrid λ=0"), ("hybrid", 0.25, "λ=0.25"), ("hybrid", 0.5, "λ=0.5"),
             ("hybrid", 0.75, "λ=0.75"), ("hybrid", 1.0, "λ=1")]
same_lambda(a, b) = (isnan(a) && isnan(b)) || a == b

panels = Plots.Plot[]
for (row, (e, sc, title)) in enumerate(PANELS)
    sub = filter(r -> r.n in (12, 40), panel(layout, e, sc))
    # a handful of layouts on the gas relative scale fall near -7 (a control on a dead
    # well); clip the axis at the 1st percentile so the rest stay readable, and say so
    lo, hi = quantile(sub.skill, 0.01), maximum(sub.skill); pad = 0.06 * (hi - lo)
    for n in (12, 40)
        p = plot(title = "$title, n = $n", ylims = (lo - pad, hi + pad), legend = false,
                 ylabel = n == 12 ? "GP skill per layout" : "",
                 xticks = (1:length(OBJ), row == length(PANELS) ? last.(OBJ) : fill("", length(OBJ))),
                 xrotation = 40, xlims = (0.4, length(OBJ) + 0.6), grid = :y)
        hline!(p, [0], color = RULE, lw = 1)
        for (k, (o, λ, _)) in enumerate(OBJ)
            v = filter(r -> r.n == n && r.objective == o && same_lambda(r.lambda, λ), sub).skill
            jit = 0.25 .* (2 .* rand(MersenneTwister(k), length(v)) .- 1)
            scatter!(p, k .+ jit, v, color = o == "random" ? RULE : C1,
                     alpha = 0.55, ms = 3, markerstrokewidth = 0)
            plot!(p, [k - 0.32, k + 0.32], fill(mean(v), 2), color = INK, lw = 2)
        end
        off = filter(r -> r.n == n && r.skill < lo - pad, sub).skill
        isempty(off) || annotate!(p, length(OBJ) + 0.5, lo - pad + 0.04 * (hi - lo),
            text(@sprintf("%d layout%s below the axis (min %.1f)", length(off), length(off) == 1 ? "" : "s", minimum(off)),
                 7, MUTED, :right))
        push!(panels, p)
    end
end
save(plot(panels..., layout = (length(PANELS), 2), size = (1000, 1180), left_margin = 4Plots.mm,
          bottom_margin = 3Plots.mm), "fig_exp1_placement")

# --- Experiment 2a: interval coverage ------------------------------------------
p = plot(title = "GP 95% interval coverage", xticks = NTICKS,
         xlims = (0, 60), ylims = (0, 1.02), yticks = 0:0.2:1,
         xlabel = "number of control wells", ylabel = "fraction of held-out wells inside interval",
         legend = false, size = (780, 470), left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)
hline!(p, [0.95], color = MUTED, lw = 1, ls = :dash)
annotate!(p, 0.5, 0.975, text("nominal 0.95", 7, MUTED, :left))
ends = Float64[]
for ((e, sc, title), c, shape) in zip(PANELS, (C1, C2, C3, C4), (:circle, :utriangle, :rect, :diamond))
    d = sort(filter(r -> r.model == "gp", panel(by_n, e, sc)), :n)
    # dotted: the 2026-09-04 version, smooth-surface sd only (per-well noise left out)
    plot!(p, d.n, d.cover95_f_mean, color = c, lw = 1, ls = :dot, alpha = 0.7)
    plot!(p, d.n, d.cover95_mean, color = c, lw = 2, marker = shape, ms = 4, markerstrokewidth = 0)
    push!(ends, last(d.cover95_mean))
end
for ((_, _, title), y) in zip(PANELS, spread(ends .- 0.04, 0.04))
    annotate!(p, 41, y, text(replace(title, ", for comparison" => ""), 7, INK, :left))
end
annotate!(p, 0.5, 0.05, text("solid: predictive sd (noise included)   dotted: noise left out (2026-09-04 version)", 7, MUTED, :left))
save(p, "fig_exp2_coverage")

# --- Experiment 2b: does posterior sd predict realised error? ------------------
panels = Plots.Plot[]
for (row, (e, sc, title)) in enumerate(PANELS), n in (12, 40)
    d = filter(r -> r.n == n && isfinite(r.post_sd) && isfinite(r.rmse), panel(layout, e, sc))
    ρ = corspearman(d.post_sd, d.rmse)
    local p = scatter(d.post_sd, d.rmse, color = C1, alpha = 0.45, ms = 3, markerstrokewidth = 0,
                legend = false, title = @sprintf("%s, n = %d  ·  ρ = %.2f", title, n, ρ),
                titlefontsize = 9,
                xlabel = "mean GP predictive sd ($(UNIT[sc]))",
                ylabel = n == 12 ? "realised RMSE ($(UNIT[sc]))" : "")
    push!(panels, p)
end
save(plot(panels..., layout = (length(PANELS), 2), size = (950, 1200), left_margin = 4Plots.mm,
          bottom_margin = 2Plots.mm), "fig_exp2_sd_vs_error")
