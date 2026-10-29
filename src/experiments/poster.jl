# The 2021 poster plate (E4): does correction make replicate wells agree?
#
# E4 is a regular 384-well plate, not an all-control plate, so there is no measured
# error surface to score against. Its design supplies a different, fair test: 24
# no-glucose and 24 top-glucose controls train the correction, and 49 middle glucose
# levels are each spread over 6-7 wells across the plate. Those middle wells are never
# used for fitting, so if correction removes position effects, wells that received the
# same glucose should agree more closely afterwards.

"""
    read_mantis_layout(path) -> NamedTuple

Glucose concentration per well from a Mantis dispense file (one row per well, with
explicit `Row` and `Col` columns), plus the control and middle-level masks:
negatives are the lowest concentration, positives the highest.
"""
function read_mantis_layout(path::AbstractString)
    rows = CSV.File(path)
    nr, nc = maximum(r.Row for r in rows), maximum(r.Col for r in rows)
    conc = fill(NaN, nr, nc)
    for r in rows
        conc[r.Row, r.Col] = Float64(r.Final_Glucose_Conc)
    end
    any(isnan, conc) && error("$(basename(path)): layout does not cover every well")
    lo, hi = extrema(conc)
    neg, pos = conc .== lo, conc .== hi
    return (conc = conc, negatives = BitMatrix(neg), positives = BitMatrix(pos),
            middle = BitMatrix(.!(neg .| pos)))
end

"""Apply a predicted error `e` on its correction's scale: subtract an offset, or remove a fraction."""
apply_correction(x, e, correction::Function) =
    correction === absolute ? x .- e : x .* (1 .- e)

"""Wells-from-the-edge index: 1 for the outer ring, rising toward the plate centre."""
ring_index(nr, nc) = [min(i, j, nr + 1 - i, nc + 1 - j) for i in 1:nr, j in 1:nc]

const POSTER_METHODS = ("uncorrected", "subtract_and_normalize", "gp_iso", "gp_ard",
                        "idw", "rowcol", "nearest")

# GP configurations, both in well units (pitch 1) like the original package. `:ard`
# reproduces PlateErrorGP before its "isotropic kernel" commit: start (column 4, row 3).
const POSTER_GP = Dict("gp_iso" => (kernel = :iso, ls = 3.0, lb_ls = 2.0, ub_ls = 6.0),
                       "gp_ard" => (kernel = :ard, ls = (4.0, 3.0), lb_ls = 2.0, ub_ls = 6.0))

"""
    two_pass(G, neg, pos, method) -> (corrected, info)

Correct a regular plate the way `gp_correction` does: negatives remove an additive
error, then positives remove a relative one from the pass-1 output. Every method uses
the same two passes and the same control wells, so they differ only in how the error
between controls is predicted. `info` records fit diagnostics for the GP methods.
"""
function two_pass(G::AbstractMatrix, neg::BitMatrix, pos::BitMatrix, method::AbstractString)
    G = Float64.(G)
    info = Dict{String,Any}()
    if method == "uncorrected"
        return G, info
    elseif method == "subtract_and_normalize"
        return subtract_and_normalize(G, PlateArray(trues(size(G)), pos, neg)), info
    elseif haskey(POSTER_GP, method)
        cfg = POSTER_GP[method]
        opts = (verbose = true, pitch = 1.0, kernel = cfg.kernel, ls = cfg.ls,
                lb_ls = cfg.lb_ls, ub_ls = cfg.ub_ls)
        c1, d1 = plateerrorGP(G, neg, absolute; opts...)
        c2, d2 = plateerrorGP(c1, pos, relative; opts...)
        for (tag, d, mask, corr, X) in (("p1", d1, neg, absolute, G), ("p2", d2, pos, relative, c1))
            gp = d["gp"]
            ell = lengthscales(gp)
            info["$(tag)_ls"] = ell
            info["$(tag)_logparams"] = GaussianProcesses.get_params(gp.kernel)
            info["$(tag)_noise_var"] = exp(2 * gp.logNoise.value)
            info["$(tag)_at_bound"] = any(l -> isapprox(l, cfg.lb_ls; rtol = 1e-3) ||
                                               isapprox(l, cfg.ub_ls; rtol = 1e-3), ell)
            y = truth_surface(X, mask, corr)[mask]
            info["$(tag)_sd_fallback"] = !(isfinite(std(y)) && std(y) > 0)
            info["$(tag)_mean_error"] = d["mean_error"]
        end
        return c2, info
    elseif haskey(baseline_methods, method)
        f = baseline_methods[method]
        c1 = apply_correction(G, f(G, neg, absolute), absolute)
        c2 = apply_correction(c1, f(c1, pos, relative), relative)
        return c2, info
    end
    throw(ArgumentError("unknown method: $method"))
end

"""Within-level and total sums of squares of `Y` over the `middle` wells, grouped by `conc`."""
function within_total_ss(Y, conc, middle)
    v, levels = Y[middle], conc[middle]
    sst = sum(abs2, v .- mean(v))
    ssw = sum(sum(abs2, v[levels .== l] .- mean(v[levels .== l])) for l in unique(levels))
    return ssw, sst
end

"""
    replicate_scores(X, raw, conc, middle, neg, pos) -> NamedTuple

How well wells that received the same glucose agree, on the middle levels only.

* `eta2` - share of the middle wells' variance explained by glucose level.
* `replicate_skill` - share of the *non-glucose* variance removed relative to the raw
  plate: `1 - (SSW/SST)_corrected / (SSW/SST)_raw`. It is unchanged by any shift or
  rescaling of the plate, so it is exactly 0 for raw data and for any affine
  normalisation, 1 if replicates agree perfectly, and negative if correction scatters them.
* `within_sd` - pooled within-level sd in control-window units, i.e. as a fraction of
  (positive median - negative median), so plates and metrics share a scale.
* `centre_excess` - mean residual (well minus its level's median) in the inner rings
  (5 or more wells from the edge) minus that in the outer ring, in window units.
"""
function replicate_scores(X, raw, conc, middle, neg, pos)
    ssw, sst = within_total_ss(X, conc, middle)
    ssw0, sst0 = within_total_ss(raw, conc, middle)
    window = median(X[pos]) - median(X[neg])
    resid = similar(X, Float64); fill!(resid, NaN)
    for l in unique(conc[middle])
        idx = middle .& (conc .== l)
        resid[idx] .= (X[idx] .- median(X[idx])) ./ window
    end
    ring = ring_index(size(X)...)
    inner, outer = middle .& (ring .>= 5), middle .& (ring .== 1)
    nlev = length(unique(conc[middle]))
    return (eta2 = 1 - ssw / sst,
            replicate_skill = 1 - (ssw / sst) / (ssw0 / sst0),
            within_sd = sqrt(ssw / (count(middle) - nlev)) / window,
            centre_excess = mean(resid[inner]) - mean(resid[outer]),
            window = window, resid = resid)
end

"""
    poster_metrics(k) -> Dict{String,Matrix{Float64}}

The per-well measurements the poster tests use: OD at read 90 (what the poster
figure plotted) and each well's maximum OD over the run, with the read it peaked at.
"""
function poster_metrics(k)
    V = Float64.(k.values)
    nr, nc, _ = size(V)
    maxod = [maximum(V[i, j, :]) for i in 1:nr, j in 1:nc]
    peak = [argmax(V[i, j, :]) for i in 1:nr, j in 1:nc]
    return Dict("od_read90" => V[:, :, 90], "max_od" => maxod), peak
end

"""
    control_subsets(neg, pos, ks, nrep; seed) -> Vector{NamedTuple}

Random subsets of the plate's own control wells, `k` of each class, `nrep` per `k`
(one when `k` is the full set). Drawn once and reused for every method, kernel and
metric, so comparisons between them are paired.
"""
function control_subsets(neg::BitMatrix, pos::BitMatrix, ks, nrep; seed = 20261009)
    nidx, pidx = findall(neg), findall(pos)
    out = NamedTuple[]
    for k in ks, rep in 1:(k >= min(length(nidx), length(pidx)) ? 1 : nrep)
        rng = Random.Xoshiro(seed + 1000k + rep)
        n_k = falses(size(neg)); p_k = falses(size(pos))
        n_k[nidx[randperm(rng, length(nidx))[1:k]]] .= true
        p_k[pidx[randperm(rng, length(pidx))[1:k]]] .= true
        push!(out, (k = k, rep = rep, negatives = BitMatrix(n_k), positives = BitMatrix(p_k)))
    end
    return out
end

"""Well name such as `"B12"` for a plate position."""
wellname(ci::CartesianIndex{2}) = string(Char('A' + ci[1] - 1), ci[2])
