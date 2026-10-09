# Figures for the three experiments, drawn from the tables analyze.jl writes.
#
# Run after experiments/analyze.jl:
#   julia --project experiments/figures.jl
#
# Writes to results/ (PNG for viewing, PDF for the manuscript):
#   fig_exp1_placement      GP skill per layout, by placement objective (n = 12, 40)
#   fig_exp2_coverage       GP 95% interval coverage vs number of controls
#   fig_exp2_sd_vs_error    posterior sd vs realised RMSE, one point per layout
#
# Posterior sd is the latent-function sd (`predict_f`), without the noise term, so
# it is not expected to match RMSE in level; Exp 2 asks whether it ranks layouts.
#   fig_exp3_skill_vs_n     skill vs number of controls, GP against simple baselines
#
# All panels use absolute correction. Relative correction is in the summary tables.

using CSV, DataFrames, Statistics, StatsBase, Plots, Printf, Random
gr()

const OUT = joinpath(@__DIR__, "..", "results")

const PLATES = [("negative", "All-negative plate"),
                ("positive", "All-positive plate"),
                ("gas",      "Gas-gradient plates")]

# Categorical slots, fixed order (validated: adjacent and 3-way all-pairs CVD-safe).
const C1, C2, C3, C4 = colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a", colorant"#eda100"
const INK, MUTED, RULE = colorant"#0b0b0b", colorant"#52514e", colorant"#8a8985"

default(fontfamily = "sans-serif", framestyle = :axes, grid = true, gridalpha = 0.12,
        gridlinewidth = 0.5, foreground_color_axis = RULE, foreground_color_border = RULE,
        tickfontcolor = MUTED, guidefontcolor = INK, titlefontcolor = INK,
        titlefontsize = 11, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8,
        legend_foreground_color = nothing, background_color = :white, dpi = 200)

const NTICKS = ([3, 6, 12, 24, 40], ["3", "6", "12", "24", "40"])

function save(p, name)
    savefig(p, joinpath(OUT, name * ".png"))
    savefig(p, joinpath(OUT, name * ".pdf"))
    @info "wrote" figure = name
end

by_n   = filter(r -> r.correction == "absolute", CSV.read(joinpath(OUT, "summary_by_n.csv"), DataFrame))
layout = filter(r -> r.correction == "absolute", CSV.read(joinpath(OUT, "per_layout_gp.csv"), DataFrame))
surf   = filter(r -> r.correction == "absolute", CSV.read(joinpath(OUT, "surfaces.csv"), DataFrame))

# Ceiling: no smooth interpolator recovers the nugget, so skill is bounded by the
# structured share of variance (median over a plate type's surfaces).
ceiling = Dict(e => median(filter(r -> r.experiment == e, surf).structured_share) for (e, _) in PLATES)

# --- Experiment 3: skill vs number of controls ---------------------------------
const MODELS = [("gp", "GP (this method)", C1, :circle),
                ("idw", "inverse-distance", C2, :utriangle),
                ("rowcol", "row + column", C3, :rect),
                ("nearest", "nearest control", C4, :diamond)]
const YLIM3 = (-0.7, 0.85)

panels = map(enumerate(PLATES)) do (i, (e, title))
    p = plot(title = title, xscale = :log10, xticks = NTICKS, xlims = (2.6, 58), ylims = YLIM3,
             xlabel = "number of control wells", ylabel = i == 1 ? "skill vs. control median" : "",
             legend = i == 1 ? :topleft : false)
    hline!(p, [0], color = RULE, lw = 1, label = "")
    hline!(p, [ceiling[e]], color = MUTED, lw = 1, ls = :dash, label = "")
    annotate!(p, 2.75, ceiling[e] + 0.035, text(@sprintf("ceiling %.2f", ceiling[e]), 7, MUTED, :left))
    noted = 0
    for (m, lab, c, shape) in MODELS
        d = sort(filter(r -> r.experiment == e && r.model == m, by_n), :n)
        plot!(p, d.n, d.skill_mean, color = c, lw = 2, marker = shape, ms = 4,
              markerstrokewidth = 0, label = lab)
        if minimum(d.skill_mean) < YLIM3[1]
            annotate!(p, 55, YLIM3[1] + 0.05 + 0.07noted,
                      text(@sprintf("%s: min %.1f, off-scale", lab, minimum(d.skill_mean)), 7, MUTED, :right))
            noted += 1
        end
    end
    p
end
save(plot(panels..., layout = (1, 3), size = (1150, 400), left_margin = 4Plots.mm,
          bottom_margin = 5Plots.mm), "fig_exp3_skill_vs_n")

# --- Experiment 1: placement objectives ----------------------------------------
const OBJ = [("random", NaN, "random"), ("LHS", NaN, "LHS"), ("minimax", NaN, "minisum*"),
             ("hybrid", 0.0, "hybrid λ=0"), ("hybrid", 0.25, "λ=0.25"), ("hybrid", 0.5, "λ=0.5"),
             ("hybrid", 0.75, "λ=0.75"), ("hybrid", 1.0, "λ=1")]
same_lambda(a, b) = (isnan(a) && isnan(b)) || a == b

panels = Plots.Plot[]
for (row, (e, title)) in enumerate(PLATES)
    sub = filter(r -> r.experiment == e && r.n in (12, 40), layout)
    lo, hi = extrema(sub.skill); pad = 0.06 * (hi - lo)
    for n in (12, 40)
        p = plot(title = "$title, n = $n", ylims = (lo - pad, hi + pad), legend = false,
                 ylabel = n == 12 ? "GP skill per layout" : "",
                 xticks = (1:length(OBJ), row == length(PLATES) ? last.(OBJ) : fill("", length(OBJ))),
                 xrotation = 40, xlims = (0.4, length(OBJ) + 0.6), grid = :y)
        hline!(p, [0], color = RULE, lw = 1)
        for (k, (o, λ, _)) in enumerate(OBJ)
            v = filter(r -> r.n == n && r.objective == o && same_lambda(r.lambda, λ), sub).skill
            jit = 0.25 .* (2 .* rand(MersenneTwister(k), length(v)) .- 1)
            scatter!(p, k .+ jit, v, color = o == "random" ? RULE : C1,
                     alpha = 0.55, ms = 3, markerstrokewidth = 0)
            plot!(p, [k - 0.32, k + 0.32], fill(mean(v), 2), color = INK, lw = 2)
        end
        push!(panels, p)
    end
end
save(plot(panels..., layout = (3, 2), size = (1000, 900), left_margin = 4Plots.mm,
          bottom_margin = 3Plots.mm), "fig_exp1_placement")

# --- Experiment 2a: interval coverage ------------------------------------------
gp = filter(r -> r.model == "gp", by_n)
p = plot(title = "GP 95% interval coverage", xscale = :log10, xticks = NTICKS,
         xlims = (2.6, 70), ylims = (0, 1.02), yticks = 0:0.2:1,
         xlabel = "number of control wells", ylabel = "fraction of held-out wells inside interval",
         legend = :bottomright, size = (620, 430), left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)
hline!(p, [0.95], color = MUTED, lw = 1, ls = :dash, label = "")
annotate!(p, 2.75, 0.975, text("nominal 0.95", 7, MUTED, :left))
for ((e, title), c, shape) in zip(PLATES, (C1, C2, C3), (:circle, :utriangle, :rect))
    d = sort(filter(r -> r.experiment == e, gp), :n)
    plot!(p, d.n, d.cover95_mean, color = c, lw = 2, marker = shape, ms = 4,
          markerstrokewidth = 0, label = title)
    annotate!(p, 44, last(d.cover95_mean), text(split(title)[1], 7, INK, :left))
end
save(p, "fig_exp2_coverage")

# --- Experiment 2b: does posterior sd predict realised error? ------------------
panels = Plots.Plot[]
for (row, (e, title)) in enumerate(PLATES), n in (12, 40)
    d = filter(r -> r.experiment == e && r.n == n && isfinite(r.post_sd) && isfinite(r.rmse), layout)
    ρ = corspearman(d.post_sd, d.rmse)
    local p = scatter(d.post_sd, d.rmse, color = C1, alpha = 0.45, ms = 3, markerstrokewidth = 0,
                legend = false, title = @sprintf("%s, n = %d  ·  ρ = %.2f", first(split(title)), n, ρ),
                xlabel = row == length(PLATES) ? "mean GP posterior sd (OD)" : "",
                ylabel = n == 12 ? "realised RMSE (OD)" : "")
    push!(panels, p)
end
save(plot(panels..., layout = (3, 2), size = (900, 900), left_margin = 4Plots.mm,
          bottom_margin = 2Plots.mm), "fig_exp2_sd_vs_error")
