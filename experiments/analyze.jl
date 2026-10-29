# Summarise the sweep into tables suitable for inspection and for figures.
#
# The unit of replication is the LAYOUT, not the layout-timepoint pair: timepoints
# within one kinetic run are autocorrelated, so a layout's score is averaged over
# surfaces first and objectives are compared across layouts.
#
# Writes to results/:
#   summary_by_n.csv          accuracy vs number of controls, per model
#   summary_by_objective.csv  layout objectives compared at each control count
#   summary_link.csv          does posterior sd predict realised error?
#   summary_by_surface.csv    per-surface GP skill against that surface's ceiling
#   per_layout_gp.csv         GP scores per layout, averaged over surfaces (for figures)

using CSV, DataFrames, Statistics, StatsBase

const OUT = joinpath(@__DIR__, "..", "results")
const COLS = [:surface, :plate, :metric, :n, :objective, :lambda, :layout, :model,
              :skill, :rmse, :mae, :bias, :spearman, :cover68, :cover95,
              :mean_sd, :max_sd, :cover95_f, :mean_sd_f, :skill_screened,
              :cover95_screened, :n_screened, :ls_fit, :noise_fit, :max_min_dist, :mean_min_dist,
              :cols_covered, :n_heldout]

corr_of(s) = endswith(s, "_abs") ? "absolute" : "relative"
nanmean(v) = (w = collect(skipmissing(v)); w = filter(isfinite, w); isempty(w) ? NaN : mean(w))
nanstd(v)  = (w = collect(skipmissing(v)); w = filter(isfinite, w); length(w) < 2 ? NaN : std(w))
nanmed(v)  = (w = collect(skipmissing(v)); w = filter(isfinite, w); isempty(w) ? NaN : median(w))

function load(exp)
    p = joinpath(OUT, "sweep_$(exp).csv")
    isfile(p) || return nothing
    df = CSV.read(p, DataFrame; select = COLS)
    df.experiment .= exp
    df.correction = corr_of.(df.surface)
    return df
end

dfs = filter(!isnothing, [load(e) for e in ("negative", "positive", "gas")])
isempty(dfs) && error("no sweep output found in $OUT")
df = vcat(dfs...)
@info "loaded" rows = nrow(df) experiments = unique(df.experiment)

# --- accuracy vs control count -------------------------------------------------
by_n = combine(groupby(df, [:experiment, :correction, :n, :model]),
               :skill => nanmean => :skill_mean,
               :skill => nanmed => :skill_median,
               :skill => nanstd => :skill_sd,
               :rmse => nanmean => :rmse_mean,
               :cover95 => nanmean => :cover95_mean,
               :cover95_f => nanmean => :cover95_f_mean,
               :cover95_screened => nanmean => :cover95_screened_mean,
               :skill_screened => nanmean => :skill_screened_mean,
               :mean_sd => nanmean => :post_sd_mean,
               :mean_sd_f => nanmean => :post_sd_f_mean,
               :ls_fit => nanmed => :ls_median,
               nrow => :n_obs)
sort!(by_n, [:experiment, :correction, :model, :n])
CSV.write(joinpath(OUT, "summary_by_n.csv"), by_n)

# --- objectives compared, layout as the unit -----------------------------------
gp = filter(r -> r.model == "gp", df)
per_layout = combine(groupby(gp, [:experiment, :correction, :n, :objective, :lambda, :layout]),
                     :skill => nanmean => :skill,
                     :rmse => nanmean => :rmse,
                     :skill_screened => nanmean => :skill_screened,
                     :mean_sd => nanmean => :post_sd,
                     :mean_sd_f => nanmean => :post_sd_f,
                     :cover95 => nanmean => :cover95,
                     :cover95_f => nanmean => :cover95_f,
                     :max_min_dist => first => :max_min_dist,
                     :cols_covered => first => :cols_covered)
CSV.write(joinpath(OUT, "per_layout_gp.csv"), per_layout)
by_obj = combine(groupby(per_layout, [:experiment, :correction, :n, :objective, :lambda]),
                 :skill => nanmean => :skill_mean,
                 :skill => nanstd => :skill_sd,
                 :rmse => nanmean => :rmse_mean,
                 :post_sd => nanmean => :post_sd_mean,
                 :cover95 => nanmean => :cover95_mean,
                 :max_min_dist => nanmean => :max_min_dist,
                 :cols_covered => nanmean => :cols_covered,
                 nrow => :n_layouts)
sort!(by_obj, [:experiment, :correction, :n, :objective, :lambda])
CSV.write(joinpath(OUT, "summary_by_objective.csv"), by_obj)

# --- does the design criterion predict realised error? -------------------------
# The manuscript claim rests on this: ranking layouts by posterior sd is only
# interesting if that ranking matches the error actually observed.
function rankcor(a, b)
    keep = isfinite.(a) .& isfinite.(b)
    sum(keep) < 4 && return NaN
    corspearman(collect(a[keep]), collect(b[keep]))
end
link = combine(groupby(per_layout, [:experiment, :correction, :n]),
               [:post_sd, :rmse] => ((a, b) -> rankcor(a, b)) => :spearman_sd_vs_rmse,
               [:post_sd_f, :rmse] => ((a, b) -> rankcor(a, b)) => :spearman_sdf_vs_rmse,
               [:max_min_dist, :rmse] => ((a, b) -> rankcor(a, b)) => :spearman_dist_vs_rmse,
               [:post_sd, :skill] => ((a, b) -> rankcor(a, b)) => :spearman_sd_vs_skill,
               nrow => :n_layouts)
sort!(link, [:experiment, :correction, :n])
CSV.write(joinpath(OUT, "summary_link.csv"), link)

# --- per surface, against its own ceiling --------------------------------------
surf = CSV.read(joinpath(OUT, "surfaces.csv"), DataFrame)
bys = combine(groupby(filter(r -> r.model == "gp" && r.n == 12, df),
                      [:experiment, :correction, :surface]),
              :skill => nanmean => :gp_skill_n12,
              :cover95 => nanmean => :cover95,
              :skill_screened => nanmean => :gp_skill_screened_n12,
              :ls_fit => nanmed => :ls_median)
bys = leftjoin(bys, select(surf, [:surface, :cv_pct, :structured_share, :rowcol_share, :morans_I]),
               on = :surface)
bys.frac_of_ceiling = bys.gp_skill_n12 ./ bys.structured_share
sort!(bys, [:experiment, :correction, :surface])
CSV.write(joinpath(OUT, "summary_by_surface.csv"), bys)

for f in ("summary_by_n.csv", "summary_by_objective.csv", "summary_link.csv", "summary_by_surface.csv")
    @info "wrote" file = f bytes = filesize(joinpath(OUT, f))
end
