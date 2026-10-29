# Experiment 3 extended to the full plate: skill vs number of controls at every
# control count from 3 (the GP's minimum) to 96.
#
# The sweep stops at n = 40 because PlateArrays places both control classes, so 2n
# must stay below 96. Here only the class being scored is placed, using uniformly
# random layouts (the sweep's "random" objective draws from the same distribution),
# so n can reach 96. Each count gets 30 layouts, shared across surfaces as in the
# sweep. The exceptions are n = 95, which gets all 96 leave-one-out layouts, and
# n = 96, which has only one layout.
#
# Every model is scored twice. The first score uses held-out wells, as the sweep
# does; at n = 96 there are none. The second uses the control wells the model was
# fit to (in-sample). Rows store sums of squares rather than ratios, so skill can
# be pooled across layouts. With one or two held-out wells, a per-layout ratio
# measures how far those few wells sit from the median more than it measures the model.
#
#   julia --project -t 64 experiments/run_full_n.jl [ns] [nlayouts]
#
# Writes results/full_n/full_n.csv.

using LabQC, CSV, DataFrames, Statistics, Random, LinearAlgebra, Dates

BLAS.set_num_threads(1)   # one fit per thread; threaded BLAS on 96x96 only oversubscribes

const RAW = joinpath(@__DIR__, "..", "raw_data")
const OUT = joinpath(@__DIR__, "..", "results", "full_n")
# the panels of fig_exp3_skill_vs_n: (experiment, correction tag)
const PANELS = [("negative", "abs"), ("positive", "rel"), ("gas", "rel"), ("gas", "abs")]
const NS = length(ARGS) >= 1 ? parse.(Int, split(ARGS[1], ",")) : collect(3:96)
const NLAYOUTS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 30
const SEED = 20261009

"""Control masks for `n` controls on a 96-well plate: random, leave-one-out at 95, the full plate at 96."""
function masks_for(n; nlayouts = NLAYOUTS)
    nw = 96
    n == nw && return [trues(8, 12)]
    n == nw - 1 && return [(m = trues(8, 12); m[k] = false; m) for k in 1:nw]
    rng = MersenneTwister(SEED + n)
    return [(m = falses(8, 12); m[randperm(rng, nw)[1:n]] .= true; m) for _ in 1:nlayouts]
end

"""Fit the GP and every baseline to one layout; return sums of squares on held-out and control wells."""
function score(s::Surface, mask::BitMatrix, pitch)
    truth = truth_surface(s.values, mask, s.correction)
    _, diag = LabQC.plateerrorGP(s.values, mask, s.correction; verbose = true, pitch = pitch,
                                 ls = 3 * pitch, lb_ls = 2 * pitch, ub_ls = 6 * pitch)
    gp = diag["gp"]
    preds = Dict{String,Matrix{Float64}}("gp" => diag["mean_error"])
    for (name, f) in baseline_methods
        preds[name] = f(s.values, mask, s.correction)
    end
    # predictive sd: held-out truth is an observed well, noise included
    sd_y = sqrt.(diag["var_y"])
    held = .!mask
    rows = NamedTuple[]
    for (model, P) in preds
        isgp = model == "gp"
        e = P .- truth
        push!(rows, (model = model, n_held = count(held),
                     sse_held = sum(abs2, e[held]), sst_held = sum(abs2, truth[held]),
                     sse_ctrl = sum(abs2, e[mask]), sst_ctrl = sum(abs2, truth[mask]),
                     in68 = isgp ? count(abs.(e[held]) .<= 1.0 .* sd_y[held]) : -1,
                     in95 = isgp ? count(abs.(e[held]) .<= 1.96 .* sd_y[held]) : -1,
                     ls_fit = isgp ? lengthscales(gp)[1] : NaN,
                     sd_fit = isgp ? signal_sd(gp) : NaN,
                     noise_fit = isgp ? exp(gp.logNoise.value) : NaN))
    end
    return rows
end

mkpath(OUT)
@info "full-n run starting" time = now() threads = Threads.nthreads() ns = (first(NS), last(NS), length(NS)) nlayouts = NLAYOUTS

surfaces = build_surfaces(RAW)
sel = [(e, s) for (e, s) in surfaces if (e, s.meta["correction"]) in PANELS]
@info "surfaces" total = length(sel)

bank = Dict(n => masks_for(n) for n in NS)
pitch = well_pitch(8, 12)
tasks = [(si, n, li) for si in eachindex(sel) for n in NS for li in eachindex(bank[n])]
@info "fits to run" tasks = length(tasks)

results = Vector{Vector{NamedTuple}}(undef, length(tasks))
failed = Threads.Atomic{Int}(0)
done = Threads.Atomic{Int}(0)
t0 = time()
Threads.@threads :dynamic for i in eachindex(tasks)
    si, n, li = tasks[i]
    e, s = sel[si]
    rows = try
        score(s, bank[n][li], pitch)
    catch err
        Threads.atomic_add!(failed, 1)
        @warn "fit failed" surface = s.id n layout = li exception = err
        NamedTuple[]
    end
    results[i] = [merge((experiment = e, correction = s.meta["correction"], surface = s.id,
                         n = n, layout = li), r) for r in rows]
    k = Threads.atomic_add!(done, 1) + 1
    k % max(1, length(tasks) ÷ 20) == 0 &&
        @info "progress" done = k of = length(tasks) elapsed_s = round(time() - t0, digits = 1)
end

df = DataFrame(reduce(vcat, results))
CSV.write(joinpath(OUT, "full_n.csv"), df)
@info "full-n run complete" rows = nrow(df) failed = failed[] seconds = round(time() - t0, digits = 1)
