# Figures for the full-n run (experiments/run_full_n.jl): Experiment 3 and the
# coverage panel of Experiment 2, with the control count running to the full plate.
#
#   julia --project experiments/figures_full_n.jl
#
# Reads results/full_n/full_n.csv, plus the sweep's summary_by_n.csv and surfaces.csv
# for comparison. Writes to results/full_n/:
#   summary_full_n.csv     skill, in-sample skill and coverage by panel, n and model
#   fig_full_n_skill       held-out skill (top) and in-sample skill (bottom) vs n
#   fig_full_n_coverage    GP 95% interval coverage vs n
#
# Held-out skill is pooled within each surface (1 - total SSE / total SS over its
# layouts) and then averaged across surfaces, so surfaces keep equal weight as in the
# sweep. Pooling matters near the full plate. At n = 95 a layout holds out a single
# well, so a per-layout ratio mostly measures how far that well sits from the median.
# The sweep's per-layout mean is kept in the summary table for comparison.

using CSV, DataFrames, Statistics, Plots, Printf
gr()

const OUT = joinpath(@__DIR__, "..", "results", "full_n")
const RES = joinpath(@__DIR__, "..", "results")

const PANELS = [("negative", "abs", "All-negative plate (additive)"),
                ("positive", "rel", "All-positive plate (relative)"),
                ("gas",      "rel", "Gas plates (relative)"),
                ("gas",      "abs", "Gas plates (additive, for comparison)")]
const LONG = Dict("abs" => "absolute", "rel" => "relative")

const C1, C2, C3, C4 = colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a", colorant"#eda100"
const INK, MUTED, RULE = colorant"#0b0b0b", colorant"#52514e", colorant"#8a8985"

default(fontfamily = "sans-serif", framestyle = :axes, grid = true, gridalpha = 0.12,
        gridlinewidth = 0.5, foreground_color_axis = RULE, foreground_color_border = RULE,
        tickfontcolor = MUTED, guidefontcolor = INK, titlefontcolor = INK,
        titlefontsize = 11, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8,
        legend_foreground_color = nothing, background_color = :white, dpi = 200)

const XT = 0:12:96

function save(p, name)
    savefig(p, joinpath(OUT, name * ".png"))
    savefig(p, joinpath(OUT, name * ".pdf"))
    @info "wrote" figure = name
end

nanmean(x) = (v = filter(!isnan, collect(skipmissing(x))); isempty(v) ? NaN : mean(v))
pooled(sse, sst) = sum(sst) > 0 ? 1 - sum(sse) / sum(sst) : NaN
# in68/in95 are -1 for the baselines, which have no interval
frac(k, n) = (sum(n) > 0 && all(>=(0), k)) ? sum(k) / sum(n) : NaN

df = CSV.read(joinpath(OUT, "full_n.csv"), DataFrame)
df.skill_layout = [r.sst_held > 0 ? 1 - r.sse_held / r.sst_held : NaN for r in eachrow(df)]

per_surface = combine(groupby(df, [:experiment, :correction, :surface, :n, :model]),
    [:sse_held, :sst_held] => pooled => :held,
    :skill_layout => nanmean => :held_layout,
    [:sse_ctrl, :sst_ctrl] => pooled => :insample,
    [:in95, :n_held] => frac => :cover95,
    [:in68, :n_held] => frac => :cover68,
    :ls_fit => nanmean => :ls_fit,
    :noise_fit => nanmean => :noise_fit)

summ = combine(groupby(per_surface, [:experiment, :correction, :n, :model]),
    :held => nanmean => :skill_pooled,
    :held_layout => nanmean => :skill_layout_mean,
    :insample => nanmean => :skill_insample,
    :cover95 => nanmean => :cover95,
    :cover68 => nanmean => :cover68,
    :ls_fit => nanmean => :ls_fit,
    :noise_fit => nanmean => :noise_fit,
    nrow => :n_surfaces)
sort!(summ, [:experiment, :correction, :model, :n])
CSV.write(joinpath(OUT, "summary_full_n.csv"), summ)
@info "wrote summary_full_n.csv" rows = nrow(summ)

# the sweep's dashed "ceiling" (median structured share), and its GP curve (n <= 40)
surf = CSV.read(joinpath(RES, "surfaces.csv"), DataFrame)
by_n = CSV.read(joinpath(RES, "summary_by_n.csv"), DataFrame)
lag1 = Dict((e, sc) => median(filter(r -> r.experiment == e && r.correction == LONG[sc], surf).structured_share)
            for (e, sc, _) in PANELS)
panel(d, e, sc) = filter(r -> r.experiment == e && r.correction == sc, d)

const MODELS = [("gp", "GP (this method)", C1),
                ("idw", "inverse-distance", C2),
                ("rowcol", "row + column", C3),
                ("nearest", "nearest control", C4)]
const YL = (-0.7, 1.05)

function skill_panel(i, e, sc, title, col; top)
    p = plot(title = top ? title : "", xlims = (0, 98), xticks = XT, ylims = YL,
             xlabel = top ? "" : "number of control wells",
             ylabel = i == 1 ? (top ? "held-out skill vs. control median" :
                                      "skill on the control wells (in-sample)") : "",
             legend = (i == 1 && top) ? :topleft : false)
    hline!(p, [0], color = RULE, lw = 1, label = "")
    hline!(p, [1], color = RULE, lw = 0.5, ls = :dot, label = "")
    vline!(p, [40], color = RULE, lw = 0.5, alpha = 0.6, label = "")
    if top
        c = lag1[(e, sc)]
        hline!(p, [c], color = MUTED, lw = 1, ls = :dash, label = "")
        annotate!(p, 1.5, c + 0.04, text(@sprintf("ceiling estimate %.2f (not a bound)", c), 7, MUTED, :left))
        annotate!(p, 41, YL[1] + 0.05, text("sweep stops at 40", 7, MUTED, :left))
        sw = sort(filter(r -> r.model == "gp", panel(by_n, e, LONG[sc])), :n)
        plot!(p, sw.n, sw.skill_mean, color = C1, lw = 1, ls = :dot, label = "GP, sweep figure (per-layout mean)")
    end
    noted = 0
    for (m, lab, c) in MODELS
        d = sort(filter(r -> r.model == m && isfinite(r[col]), panel(summ, e, sc)), :n)
        plot!(p, d.n, d[!, col], color = c, lw = 2, label = lab)
        lo = minimum(d[!, col])
        if lo < YL[1]
            annotate!(p, 96, YL[1] + 0.12 + 0.07noted,
                      text(@sprintf("%s: min %.1f, off-scale", lab, lo), 7, MUTED, :right))
            noted += 1
        end
    end
    # nearest and IDW are exact interpolators, so they sit on 1 and overlap
    top || annotate!(p, 95, -0.6, text("nearest and inverse-distance overlap at 1:\nthey replay the control values exactly",
                                       7, MUTED, :right))
    p
end

tops = [skill_panel(i, e, sc, t, :skill_pooled; top = true) for (i, (e, sc, t)) in enumerate(PANELS)]
bots = [skill_panel(i, e, sc, t, :skill_insample; top = false) for (i, (e, sc, t)) in enumerate(PANELS)]
save(plot(tops..., bots..., layout = (2, 4), size = (1600, 860), left_margin = 5Plots.mm,
          bottom_margin = 5Plots.mm), "fig_full_n_skill")

# --- coverage ------------------------------------------------------------------
p = plot(title = "GP 95% interval coverage, held-out wells", xlims = (0, 98), xticks = XT,
         ylims = (0, 1.02), yticks = 0:0.2:1, xlabel = "number of control wells",
         ylabel = "fraction of held-out wells inside interval", legend = :bottomright,
         size = (820, 480), left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)
hline!(p, [0.95], color = MUTED, lw = 1, ls = :dash, label = "")
annotate!(p, 1.5, 0.975, text("nominal 0.95", 7, MUTED, :left))
for ((e, sc, title), c) in zip(PANELS, (C1, C2, C3, C4))
    d = sort(filter(r -> r.model == "gp" && isfinite(r.cover95), panel(summ, e, sc)), :n)
    plot!(p, d.n, d.cover95, color = c, lw = 2, label = replace(title, ", for comparison" => ""))
end
save(p, "fig_full_n_coverage")
