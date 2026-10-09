# Driver for the LabQC placement / control-count sweep.
#
# Writes to results/:
#   surfaces.csv  - one row per measured surface, with its variance decomposition.
#                   The structured share is the ceiling on skill for that surface.
#   layouts.csv   - one row per generated layout, with geometry and objective scores
#   sweep_<experiment>.csv - one row per (surface, n, objective, layout, model)
#
# Resumable: rerunning skips cells already present in the output CSVs.

using LabQC, PlateArrays, CSV, DataFrames, Statistics, Dates

const RAW = joinpath(@__DIR__, "..", "raw_data")
const OUT = joinpath(@__DIR__, "..", "results")
# both classes are placed, so 2n must leave at least one free well on a 96-well plate
const NS = [3, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40]
const NLAYOUTS = 30
const THIN = 4

mkpath(OUT)
@info "sweep starting" time = now() ns = NS nlayouts = NLAYOUTS thin = THIN

# ---------------------------------------------------------------- surfaces
surfaces = build_surfaces(RAW; thin = THIN)
@info "surfaces built" total = length(surfaces)

srows = NamedTuple[]
for (exp, s) in surfaces
    vd = variance_decomposition(s.values)
    push!(srows, (experiment = exp, surface = s.id, plate = s.meta["plate"],
                  metric = s.meta["metric"], read = s.meta["read"], time_h = s.meta["time_h"],
                  correction = string(s.correction), class = string(s.class),
                  mean = mean(s.values), sd = std(s.values),
                  cv_pct = 100 * std(s.values) / mean(s.values),
                  min = minimum(s.values), max = maximum(s.values),
                  total_var = vd.total, nugget = vd.nugget, nugget_share = vd.nugget_share,
                  structured_share = vd.structured_share, row_share = vd.row_share,
                  col_share = vd.col_share, rowcol_share = vd.rowcol_share,
                  morans_I = vd.morans_I, morans_I_null = vd.morans_I_null))
end
CSV.write(joinpath(OUT, "surfaces.csv"), srows)
@info "wrote surfaces.csv" n = length(srows)

# ---------------------------------------------------------------- layouts
objspecs = default_objectives()
@info "building layout bank" n_specs = length(objspecs) n_counts = length(NS)
t0 = time()
lb = layout_bank(trues(8, 12), NS, objspecs; nlayouts = NLAYOUTS,
                 cache = joinpath(OUT, "layout_bank.jls"))
@info "layout bank built" seconds = round(time() - t0, digits = 1)

lrows = NamedTuple[]
for ((n, obj, lam), pas) in lb.bank, (li, pa) in enumerate(pas)
    for (cls, mask) in (("negatives", pa.negatives), ("positives", pa.positives))
        g = layout_geometry(mask, 8, 12)
        sc = layout_scores(pa, lb.bounds[n]; lambda = isnan(lam) ? 0.5 : lam)
        push!(lrows, (n = n, objective = obj, lambda = lam, layout = li, class = cls,
                      wells = join(sort([string("ABCDEFGH"[c[1]], c[2]) for c in findall(mask)]), " "),
                      max_min_dist = g.max_min_dist, mean_min_dist = g.mean_min_dist,
                      rows_covered = g.rows_covered, cols_covered = g.cols_covered,
                      max_per_row = g.max_per_row, max_per_col = g.max_per_col,
                      minisum = sc.minisum, lhs_score = sc.lhs, hybrid_score = sc.hybrid))
    end
end
CSV.write(joinpath(OUT, "layouts.csv"), lrows)
@info "wrote layouts.csv" n = length(lrows)

# ---------------------------------------------------------------- sweep
for exp in EXPERIMENTS
    ss = [s for (e, s) in surfaces if e == exp]
    isempty(ss) && continue
    path = joinpath(OUT, "sweep_$(exp).csv")
    @info "running experiment" experiment = exp surfaces = length(ss) out = basename(path)
    t = time()
    n = run_experiment(exp, ss, NS, objspecs, lb, path; nlayouts = NLAYOUTS, verbose = true)
    @info "experiment done" experiment = exp cells = n seconds = round(time() - t, digits = 1)
end

@info "sweep complete" time = now()
for f in sort(readdir(OUT))
    @info "  output" file = f bytes = filesize(joinpath(OUT, f))
end
