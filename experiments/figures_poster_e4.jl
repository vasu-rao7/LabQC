# Figures for the poster plate (E4), drawn from what run_poster_e4.jl writes.
#
#   julia --project experiments/figures_poster_e4.jl
#
# Writes to results/poster_e4/ (PNG for viewing, PDF for the manuscript):
#   fig_e4_scatter           glucose vs OD before and after correction (the poster figure)
#   fig_e4_residual_maps     each well minus its glucose level's median, on the plate
#   fig_e4_rings             that residual averaged by distance from the plate edge
#   fig_e4_controls_needed   replicate skill vs number of controls per class

using CSV, DataFrames, Statistics, Plots, Printf
gr()

const OUT = joinpath(@__DIR__, "..", "results", "poster_e4")

# Categorical slots in fixed order (validated adjacent-pair CVD-safe for 5 lines; slots
# 3-5 are below 3:1 contrast, so every line is also direct-labelled and has its own marker).
const C1, C2, C3, C4, C5 = colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a",
                           colorant"#eda100", colorant"#e87ba4"
const INK, MUTED, RULE = colorant"#0b0b0b", colorant"#52514e", colorant"#8a8985"
# diverging: blue (reads low) <-> red (reads high), neutral gray midpoint
const DIVERGING = cgrad([colorant"#184f95", colorant"#3987e5", colorant"#9ec5f4", colorant"#f0efec",
                         colorant"#f4b9b8", colorant"#e34948", colorant"#a52a2a"])

default(fontfamily = "sans-serif", framestyle = :axes, grid = true, gridalpha = 0.12,
        gridlinewidth = 0.5, foreground_color_axis = RULE, foreground_color_border = RULE,
        tickfontcolor = MUTED, guidefontcolor = INK, titlefontcolor = INK,
        titlefontsize = 10, guidefontsize = 9, tickfontsize = 8, legendfontsize = 8,
        legend_foreground_color = nothing, background_color = :white, dpi = 200)

function save(p, name)
    savefig(p, joinpath(OUT, name * ".png"))
    savefig(p, joinpath(OUT, name * ".pdf"))
    @info "wrote" figure = name
end

wells = CSV.read(joinpath(OUT, "wells_long.csv"), DataFrame)
ba    = CSV.read(joinpath(OUT, "before_after.csv"), DataFrame)
rings = CSV.read(joinpath(OUT, "rings.csv"), DataFrame)
cn    = filter(r -> r.status == "ok", CSV.read(joinpath(OUT, "controls_needed.csv"), DataFrame))

const METRICS = [("od_read90", "OD at read 90 (15 h)"), ("max_od", "max OD per well")]
const LABEL = Dict("uncorrected" => "uncorrected", "gp_iso" => "GP, one lengthscale",
                   "gp_ard" => "GP, poster kernel", "idw" => "inverse-distance",
                   "rowcol" => "row + column", "nearest" => "nearest control",
                   "subtract_and_normalize" => "subtract_and_normalize")
skill_of(m, meth) = only(filter(r -> r.metric == m && r.method == meth, ba)).replicate_skill

# Direct labels at line ends, nudged apart so lines that finish close together stay legible.
function spread(ys, gap)
    order = sortperm(ys)
    out = float.(ys[order])
    for i in 2:length(out)
        out[i] = max(out[i], out[i-1] + gap)
    end
    out .-= (sum(out) - sum(ys)) / length(out)       # keep the group centred on the data
    res = similar(out); res[order] = out
    return res
end

# --- 1. the poster figure, before and after --------------------------------------
const SCATTER_METHODS = ["uncorrected", "gp_iso", "gp_ard"]
panels = Plots.Plot[]
for (row, (m, mlabel)) in enumerate(METRICS), (col, meth) in enumerate(SCATTER_METHODS)
    d = filter(r -> r.metric == m && r.method == meth, wells)
    mid = filter(r -> r.class == "middle", d)
    inner, outer = filter(r -> r.ring >= 5, mid), filter(r -> r.ring < 5, mid)
    ctl = filter(r -> r.class != "middle", d)
    ttl = meth == "uncorrected" ? "uncorrected" :
          @sprintf("%s  ·  replicate skill %.2f", LABEL[meth], skill_of(m, meth))
    p = plot(title = ttl, legend = row == 1 && col == 1 ? :topleft : false,
             xlabel = row == 2 ? "glucose (g/L)" : "", ylabel = col == 1 ? mlabel : "",
             ylims = (0, 0.8))
    scatter!(p, 1000 .* outer.conc, outer.corrected, color = RULE, alpha = 0.6, ms = 3,
             markerstrokewidth = 0, label = "middle level, outer wells")
    scatter!(p, 1000 .* inner.conc, inner.corrected, color = C1, alpha = 0.8, ms = 3,
             markerstrokewidth = 0, label = "middle level, ≥5 wells from edge")
    scatter!(p, 1000 .* ctl.conc, ctl.corrected, color = :white, markerstrokecolor = INK,
             markerstrokewidth = 1, ms = 3.5, marker = :diamond, label = "controls (fit)")
    push!(panels, p)
end
save(plot(panels..., layout = (2, 3), size = (1200, 720), left_margin = 4Plots.mm,
          bottom_margin = 4Plots.mm), "fig_e4_scatter")

# --- 2. residual maps -----------------------------------------------------------
const MAP_METHODS = ["uncorrected", "gp_iso", "gp_ard", "idw"]
lim = maximum(abs, skipmissing(filter(r -> r.class == "middle" && r.method in MAP_METHODS, wells).resid))
lim = ceil(lim * 10) / 10
panels = Plots.Plot[]
for (row, (m, mlabel)) in enumerate(METRICS), meth in MAP_METHODS
    d = filter(r -> r.metric == m && r.method == meth, wells)
    M = fill(NaN, 16, 24)
    for r in eachrow(d)
        r.class == "middle" && (M[r.row, r.col] = r.resid)
    end
    p = heatmap(1:24, 1:16, M, yflip = true, clims = (-lim, lim), color = DIVERGING,
                title = "$(LABEL[meth]), $(m == "od_read90" ? "read 90" : "max OD")",
                yticks = (1:16, string.(collect("ABCDEFGHIJKLMNOP"))), xticks = [1, 6, 12, 18, 24],
                aspect_ratio = :equal, grid = false, colorbar = false,
                xlims = (0.5, 24.5), ylims = (0.5, 16.5), framestyle = :box)
    push!(panels, p)
end
# one shared colour bar in its own narrow column, so every map is drawn the same size
cb = scatter([0.0], [0.0], zcolor = [0.0], clims = (-lim, lim), color = DIVERGING, ms = 0,
             markerstrokewidth = 0, framestyle = :none, legend = false, colorbar = true,
             colorbar_title = "well − its level's median (control-window units);  blank = control well")
save(plot(panels..., cb, layout = @layout([grid(2, 4) a{0.07w}]), size = (1600, 560),
          left_margin = 2Plots.mm, right_margin = 6Plots.mm), "fig_e4_residual_maps")

# --- 3. residual by distance from the edge --------------------------------------
const RING_METHODS = [("uncorrected", RULE, :circle), ("gp_iso", C1, :circle),
                      ("gp_ard", C5, :utriangle), ("idw", C2, :rect)]
panels = Plots.Plot[]
for (i, (m, mlabel)) in enumerate(METRICS)
    p = plot(title = mlabel, xlabel = "wells from the plate edge (1 = outer ring)",
             ylabel = i == 1 ? "mean of (well − its level's median)\n(control-window units)" : "",
             legend = false, xlims = (0.5, 9.8), xticks = 1:8)
    hline!(p, [0], color = RULE, lw = 1)
    ends = Float64[]
    for (meth, c, shape) in RING_METHODS
        d = sort(filter(r -> r.metric == m && r.method == meth, rings), :ring)
        plot!(p, d.ring, d.mean_resid, color = c, lw = 2, marker = shape, ms = 4, markerstrokewidth = 0)
        push!(ends, last(d.mean_resid))
    end
    for ((meth, _, _), y) in zip(RING_METHODS, spread(ends, 0.022))
        annotate!(p, 8.15, y, text(LABEL[meth], 7, INK, :left))
    end
    push!(panels, p)
end
save(plot(panels..., layout = (1, 2), size = (1150, 420), left_margin = 6Plots.mm,
          bottom_margin = 5Plots.mm, right_margin = 12Plots.mm), "fig_e4_rings")

# --- 4. controls needed ---------------------------------------------------------
const CN_STYLE = [("gp_iso", C1, :circle), ("gp_ard", C5, :utriangle), ("idw", C2, :rect),
                  ("rowcol", C3, :diamond), ("nearest", C4, :star5)]
panels = Plots.Plot[]
for (i, (m, mlabel)) in enumerate(METRICS)
    p = plot(title = mlabel, xlabel = "controls per class (of 24)",
             ylabel = i == 1 ? "replicate skill\n(median and 10–90% of 50 random subsets)" : "",
             legend = false, xticks = [3, 4, 6, 8, 12, 16, 20, 24], xlims = (2, 31),
             ylims = (-0.85, 0.75))
    hline!(p, [0], color = RULE, lw = 1)
    full = only(filter(r -> r.metric == m && r.method == "gp_iso" && r.k == 24, cn)).replicate_skill
    hline!(p, [0.9full], color = MUTED, lw = 1, ls = :dash)
    annotate!(p, 2.2, 0.9full + 0.04, text("90% of the one-lengthscale GP with all 24", 7, MUTED, :left))
    ends = Float64[]
    for (meth, c, shape) in CN_STYLE
        g = combine(groupby(filter(r -> r.metric == m && r.method == meth, cn), :k),
                    :replicate_skill => median => :med,
                    :replicate_skill => (x -> quantile(x, 0.1)) => :lo,
                    :replicate_skill => (x -> quantile(x, 0.9)) => :hi)
        sort!(g, :k)
        plot!(p, g.k, g.med, ribbon = (g.med .- g.lo, g.hi .- g.med), fillalpha = 0.12,
              color = c, lw = 2, marker = shape, ms = 4, markerstrokewidth = 0)
        push!(ends, last(g.med))
    end
    for ((meth, _, _), y) in zip(CN_STYLE, spread(ends, 0.06))
        annotate!(p, 24.6, y, text(LABEL[meth], 7, INK, :left))
    end
    push!(panels, p)
end
save(plot(panels..., layout = (1, 2), size = (1250, 460), left_margin = 7Plots.mm,
          bottom_margin = 5Plots.mm, right_margin = 4Plots.mm), "fig_e4_controls_needed")
