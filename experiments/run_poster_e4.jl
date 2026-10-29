# Poster plate E4: reproduce the 2021 correction, then test whether correction makes
# replicate wells agree (fair before/after) and how many controls that takes.
#
#   julia --project experiments/run_poster_e4.jl
#
# Writes to results/poster_e4/:
#   reproduction.csv     the 2021 notebook's printed outputs vs this code (ARD, read 90)
#   before_after.csv     replicate agreement per measurement x method, all 24+24 controls
#   rings.csv            mean residual by distance from the plate edge
#   wells_long.csv       every well's raw and corrected value (for figures)
#   controls_needed.csv  replicate agreement for random subsets of k controls per class
#   control_subsets.csv  the wells in each subset
#   peak_reads.csv       the read at which each well reached its max OD

using LabQC, CSV, DataFrames, Statistics, Dates, Random
using LabQC: read_mantis_layout, two_pass, replicate_scores, poster_metrics, control_subsets,
             ring_index, wellname, POSTER_METHODS

const RAW = joinpath(@__DIR__, "..", "raw_data", "poster_e4")
const OUT = joinpath(@__DIR__, "..", "results", "poster_e4")
const KS = [3, 4, 6, 8, 12, 16, 20, 24]
const NREP = 50
mkpath(OUT)

write_run_info(OUT)    # the commit, uncommitted diff and untracked source behind these results

# ---------------------------------------------------------------- data
k = read_kinetic(joinpath(RAW, "E4_211118_CDM_Glucose_V2_cleaned.xlsx"))
lay = read_mantis_layout(joinpath(RAW, "E4_layout_SMU_Glucose_Mantis_Input2_111621.csv"))
@assert size(k.values)[1:2] == size(lay.conc) == (16, 24)
@assert count(lay.negatives) == 24 && count(lay.positives) == 24 && count(lay.middle) == 336
reps = [count(lay.middle .& (lay.conc .== l)) for l in unique(lay.conc[lay.middle])]
@assert length(reps) == 49 && all(in((6, 7)), reps)
metrics, peak = poster_metrics(k)
@info "E4 loaded" reads = k.nreads hours = k.times_h[end] levels = length(reps)

CSV.write(joinpath(OUT, "peak_reads.csv"),
          DataFrame(well = [wellname(ci) for ci in CartesianIndices(peak)][:],
                    conc = lay.conc[:], peak_read = peak[:],
                    peak_hours = k.times_h[peak[:]]))

# ---------------------------------------------------------------- reproduction
# Values printed by PlateErrorGP/GP_Notebook.ipynb at read 90: cell 56 (pass 1,
# negatives) and cell 62 (pass 2, positives), both with the original ARD kernel.
# The tolerance is set by the size of the correction itself: the raw plate already
# sits within 0.01 OD of the corrected values, so only a tight match is evidence.
const P1_ERROR = Dict("A1" => -0.00017611731471516012, "A2" => -0.0013302024899062466,
                      "A23" => -0.010698498866947292, "A24" => -0.010036779526511205,
                      "B1" => 2.3355593129579383e-5, "B24" => -0.010233480102097486,
                      "O1" => -0.003160486891780497, "O24" => -0.015477689912691286,
                      "P1" => -0.004375880900009788, "P24" => -0.010694503026647761)
const P2_VALUE = Dict("A1" => 0.15733070657529374, "A2" => 0.13543092469122742,
                      "A23" => 0.35072850378062675, "A24" => 0.43726334378147796,
                      "B1" => 0.4086355632776726, "B24" => 0.46875111483250104,
                      "O1" => 0.41599848805862094, "O24" => 0.37300362074766885,
                      "P2" => 0.0971029426287524, "P24" => 0.32016871294467675)
const P1_LOGPARAMS = [1.7917594679301678, 1.249399570246029, -2.8288877531691643]
const P2_LOGPARAMS = [1.7917594692213736, 0.693147180614958, -2.927836001984723]
const TOL_OD = 5e-4

G90 = metrics["od_read90"]
corr, info = two_pass(G90, lay.negatives, lay.positives, "gp_ard")
at(w) = (Int(w[1]) - Int('A') + 1, parse(Int, w[2:end]))
rrows = NamedTuple[]
for (w, v) in sort(collect(P1_ERROR))
    got = info["p1_mean_error"][at(w)...]
    push!(rrows, (check = "pass1 predicted error", well = w, expected = v, got = got,
                  diff = got - v, tol = TOL_OD, pass = abs(got - v) <= TOL_OD))
end
for (w, v) in sort(collect(P2_VALUE))
    got = corr[at(w)...]
    push!(rrows, (check = "pass2 corrected OD", well = w, expected = v, got = got,
                  diff = got - v, tol = TOL_OD, pass = abs(got - v) <= TOL_OD))
end
for (tag, want) in (("p1", P1_LOGPARAMS), ("p2", P2_LOGPARAMS))
    got = info["$(tag)_logparams"]
    for (name, a, b) in zip(("log ls column", "log ls row", "log signal sd"), want, got)
        push!(rrows, (check = "$tag $name", well = "", expected = a, got = b,
                      diff = b - a, tol = 0.05, pass = abs(b - a) <= 0.05))
    end
end
CSV.write(joinpath(OUT, "reproduction.csv"), rrows)
npass = count(r -> r.pass, rrows)
@info "poster reproduction" passed = npass of = length(rrows)

# ---------------------------------------------------------------- before / after
ba, rings, wells = NamedTuple[], NamedTuple[], NamedTuple[]
ringidx = ring_index(16, 24)
for (mname, G) in sort(collect(metrics); by = first), method in POSTER_METHODS
    X, inf = two_pass(G, lay.negatives, lay.positives, method)
    s = replicate_scores(X, G, lay.conc, lay.middle, lay.negatives, lay.positives)
    push!(ba, (metric = mname, method = method, replicate_skill = s.replicate_skill,
               eta2 = s.eta2, within_sd = s.within_sd, centre_excess = s.centre_excess,
               window = s.window,
               p1_ls = join(round.(get(inf, "p1_ls", Float64[]), digits = 3), "/"),
               p2_ls = join(round.(get(inf, "p2_ls", Float64[]), digits = 3), "/"),
               p1_at_bound = get(inf, "p1_at_bound", missing),
               p2_at_bound = get(inf, "p2_at_bound", missing)))
    for r in 1:maximum(ringidx)
        sel = lay.middle .& (ringidx .== r)
        any(sel) || continue
        push!(rings, (metric = mname, method = method, ring = r,
                      mean_resid = mean(s.resid[sel]), n = count(sel)))
    end
    for ci in CartesianIndices(G)
        cls = lay.negatives[ci] ? "negative" : lay.positives[ci] ? "positive" : "middle"
        push!(wells, (metric = mname, method = method, well = wellname(ci), row = ci[1],
                      col = ci[2], ring = ringidx[ci], conc = lay.conc[ci], class = cls,
                      raw = G[ci], corrected = X[ci], resid = s.resid[ci]))
    end
end
CSV.write(joinpath(OUT, "before_after.csv"), ba)
CSV.write(joinpath(OUT, "rings.csv"), rings)
CSV.write(joinpath(OUT, "wells_long.csv"), wells)
@info "before/after done" rows = length(ba)

# ---------------------------------------------------------------- controls needed
subsets = control_subsets(lay.negatives, lay.positives, KS, NREP)
CSV.write(joinpath(OUT, "control_subsets.csv"),
          [(k = s.k, rep = s.rep,
            negatives = join(wellname.(findall(s.negatives)), " "),
            positives = join(wellname.(findall(s.positives)), " ")) for s in subsets])

const CN_METHODS = ("gp_iso", "gp_ard", "idw", "rowcol", "nearest")
cn = NamedTuple[]
t0 = time()
for s in subsets, (mname, G) in sort(collect(metrics); by = first), method in CN_METHODS
    # ARD has four hyperparameters; three points cannot identify them
    if method == "gp_ard" && s.k < 4
        push!(cn, (k = s.k, rep = s.rep, metric = mname, method = method, status = "skipped",
                   error = "k < 4", replicate_skill = NaN, eta2 = NaN, within_sd = NaN,
                   centre_excess = NaN, p1_ls = "", p2_ls = "", at_bound = missing,
                   sd_fallback = missing))
        continue
    end
    try
        X, inf = two_pass(G, s.negatives, s.positives, method)
        sc = replicate_scores(X, G, lay.conc, lay.middle, lay.negatives, lay.positives)
        push!(cn, (k = s.k, rep = s.rep, metric = mname, method = method, status = "ok",
                   error = "", replicate_skill = sc.replicate_skill, eta2 = sc.eta2,
                   within_sd = sc.within_sd, centre_excess = sc.centre_excess,
                   p1_ls = join(round.(get(inf, "p1_ls", Float64[]), digits = 3), "/"),
                   p2_ls = join(round.(get(inf, "p2_ls", Float64[]), digits = 3), "/"),
                   at_bound = haskey(inf, "p1_at_bound") ? (inf["p1_at_bound"] || inf["p2_at_bound"]) : missing,
                   sd_fallback = haskey(inf, "p1_sd_fallback") ? (inf["p1_sd_fallback"] || inf["p2_sd_fallback"]) : missing))
    catch err
        push!(cn, (k = s.k, rep = s.rep, metric = mname, method = method, status = "failed",
                   error = first(sprint(showerror, err), 200), replicate_skill = NaN, eta2 = NaN,
                   within_sd = NaN, centre_excess = NaN, p1_ls = "", p2_ls = "",
                   at_bound = missing, sd_fallback = missing))
    end
end
CSV.write(joinpath(OUT, "controls_needed.csv"), cn)
@info "controls-needed done" fits = length(cn) failed = count(r -> r.status == "failed", cn) seconds = round(time() - t0)
