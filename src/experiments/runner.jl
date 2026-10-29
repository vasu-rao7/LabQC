# Sweep runner.
#
# Three plate types are run as three separate experiments, each writing its own CSV,
# so they can be examined alone or compared. Everything streams to disk as it is
# produced: a long run survives interruption, and rerunning skips cells already on
# disk rather than starting over.

const EXPERIMENTS = ("negative", "positive", "gas")

"""
    build_surfaces(rawdir; thin=4, full_timepoints=false) -> Vector{Tuple{String,Surface}}

Every measured surface to sweep, tagged with the experiment it belongs to.

* `negative` - all-negative plate 1, one surface per read (`absolute` error)
* `positive` - all-positive plate: one surface per read, plus OD at 24 h and max OD
* `gas`      - the five distinct OD600 grids from the shaker/gas series

`thin` keeps every nth kinetic read. The derived positive-plate metrics and the gas
grids are always kept. Reads inside the start-up transient are excluded: the plate is
still equilibrating there and the surface is not the error the model targets.
"""
function build_surfaces(rawdir::AbstractString; thin::Int = 4,
                        corrections = (absolute, relative))
    base = Tuple{String,Surface}[]
    out = Tuple{String,Surface}[]

    k = read_kinetic(joinpath(rawdir, "092424_neg_control_plate1.xlsx"))
    tr = last(transient_reads(k))
    for t in (tr+1):thin:k.nreads
        G = Float64.(k.values[:, :, t])
        push!(base, ("negative", Surface("neg1_t$(lpad(t,3,'0'))", G, absolute, :negatives,
              Dict{String,Any}("plate" => "neg_plate1", "metric" => "od",
                               "read" => t, "time_h" => k.times_h[t]))))
    end

    p = read_kinetic(joinpath(rawdir, "2026_08_21_VR_PosControl.xlsx"))
    ptr = last(transient_reads(p))
    for t in (ptr+1):thin:p.nreads
        G = Float64.(p.values[:, :, t])
        push!(base, ("positive", Surface("pos_t$(lpad(t,3,'0'))", G, relative, :positives,
              Dict{String,Any}("plate" => "pos_plate", "metric" => "od",
                               "read" => t, "time_h" => p.times_h[t]))))
    end
    maxod = [maximum(skipmissing(p.values[i, j, :])) for i in axes(p.values, 1), j in axes(p.values, 2)]
    push!(base, ("positive", Surface("pos_maxOD", maxod, relative, :positives,
          Dict{String,Any}("plate" => "pos_plate", "metric" => "maxOD",
                           "read" => 0, "time_h" => NaN))))
    push!(base, ("positive", Surface("pos_od24h", Float64.(p.values[:, :, end]), relative, :positives,
          Dict{String,Any}("plate" => "pos_plate", "metric" => "od24h",
                           "read" => p.nreads, "time_h" => p.times_h[end]))))

    seen = String[]
    for f in sort(filter(x -> endswith(x, ".xlsx"), readdir(joinpath(rawdir, "2026-06-04"))))
        for g in read_endpoint_grids(joinpath(rawdir, "2026-06-04", f))
            get(g.meta, "Read", "") == "OD600" || continue   # blue dye: wrong protocol, excluded
            g.hash in seen && continue                        # duplicate export of the same grid
            push!(seen, g.hash)
            push!(base, ("gas", Surface("gas_$(g.hash)", Float64.(g.values), relative, :positives,
                  Dict{String,Any}("plate" => g.source, "metric" => "od_endpoint",
                                   "read" => 0, "time_h" => NaN, "sheet" => g.sheet))))
        end
    end
    # Every surface is swept under both correction scales so the three experiments
    # are directly comparable. `relative` is unbounded as a well approaches zero,
    # which matters on the gas plates where the far columns barely grew; running
    # `absolute` alongside makes that visible instead of hiding it in an RMSE.
    for (exp, s) in base, corr in corrections
        tag = corr === absolute ? "abs" : "rel"
        push!(out, (exp, Surface(string(s.id, "_", tag), s.values, corr, s.class,
                                 merge(s.meta, Dict{String,Any}("correction" => tag)))))
    end
    return out
end

"""Objective specifications swept: the three named objectives, a random null, and hybrid at several weights."""
default_objectives(lambdas = [0.0, 0.25, 0.5, 0.75, 1.0]) =
    vcat([("random", NaN), ("minimax", NaN), ("LHS", NaN)],
         [("hybrid", l) for l in lambdas])

"""
    layout_bank(wells, ns, objspecs; nlayouts, seed) -> Dict

Pre-generate every layout once and reuse it across all surfaces, so the same control
positions are compared on every plate. Keyed by `(n, objective, lambda)`.
"""
function layout_bank(wells::BitMatrix, ns, objspecs; nlayouts = 30, seed = 20260904,
                     cache::Union{Nothing,String} = nothing)
    if cache !== nothing && isfile(cache)
        @info "layout bank loaded from cache" file = basename(cache)
        lb = deserialize(cache)
        # the cache is trusted on file existence alone, so check it was built for this run
        want = Set((n, obj, lam) for n in ns if 2n < sum(wells) for (obj, lam) in objspecs)
        Set(keys(lb.bank)) == want ||
            error("layout cache $(basename(cache)) was built for different control counts or objectives")
        all(v -> length(v) == nlayouts, values(lb.bank)) ||
            error("layout cache $(basename(cache)) does not hold $nlayouts layouts per condition")
        return lb
    end
    bounds = Dict{Int,Any}()
    bank = Dict{Tuple{Int,String,Float64},Vector{PlateArray}}()
    nwells = sum(wells)
    for n in ns
        # `exchange` swaps a control into a free well, so a layout with both classes
        # placed needs at least one well spare. PlateArrays accepts P+N == nwells at
        # its occupancy check and then fails inside the solver, so screen it here.
        if 2n >= nwells
            @warn "skipping control count: no free wells for the exchange step" n nwells
            continue
        end
        bounds[n] = fixed_hybrid_bounds(wells, n, n; reps = 4, seed = seed + n)
        for (obj, lam) in objspecs
            bank[(n, obj, lam)] = generate_layouts(wells, n, n; objective = obj,
                                                   lambda = isnan(lam) ? 0.5 : lam,
                                                   nlayouts = nlayouts, seed = seed + n)
        end
    end
    result = (bank = bank, bounds = bounds)
    cache === nothing || (mkpath(dirname(cache)); serialize(cache, result);
                          @info "layout bank cached" file = basename(cache))
    return result
end

"""
    write_run_info(outdir; repo) -> String

Record exactly which code produced the results in `outdir`: the commit, the diff
against it, and the full text of source files git does not track yet. Each start
appends to `RUN_INFO.txt` and writes its own timestamped patch, so a resumed run
keeps the record of the code that wrote its earlier rows.
"""
function write_run_info(outdir::AbstractString; repo = joinpath(@__DIR__, "..", ".."))
    stamp = Dates.format(now(), "yyyymmddTHHMMSS")
    sha = readchomp(`git -C $repo rev-parse HEAD`)
    paths = ["src", "experiments", "test", "Project.toml", "Manifest.toml"]
    io = IOBuffer()
    write(io, read(`git -C $repo diff HEAD -- $paths`, String))
    untracked = split(readchomp(`git -C $repo ls-files --others --exclude-standard -- $paths`), '\n'; keepempty = false)
    for f in untracked
        println(io, "\n=== untracked file: ", f, " ===")
        write(io, read(joinpath(repo, f), String))
    end
    patch = String(take!(io))
    name = isempty(patch) ? "none" : "RUN_DIFF_$stamp.patch"
    isempty(patch) || write(joinpath(outdir, name), patch)
    open(joinpath(outdir, "RUN_INFO.txt"), "a") do f
        println(f, "started: $(now())  commit: $sha  uncommitted changes and untracked files: $name")
    end
    return name
end

"""Write `rows` to `path`, creating the header on first use."""
function append_rows(path::AbstractString, rows::Vector{<:NamedTuple})
    isempty(rows) && return
    mkpath(dirname(path))
    if isfile(path)
        # appending under a different header would silently misalign every column
        have = String.(propertynames(CSV.File(path; limit = 1)))
        want = String.(collect(keys(first(rows))))
        have == want || error("$(basename(path)) has columns $have, new rows have $want; " *
                              "archive the old file before rerunning with a changed schema")
    end
    CSV.write(path, rows; append = isfile(path), writeheader = !isfile(path))
end

"""Keys `(surface, n, objective, lambda, layout)` already present in `path`."""
function completed_keys(path::AbstractString)
    done = Set{Tuple{String,Int,String,Float64,Int}}()
    isfile(path) || return done
    for r in CSV.File(path)
        push!(done, (String(r.surface), r.n, String(r.objective),
                     Float64(coalesce(r.lambda, NaN)), r.layout))
    end
    return done
end

"""
    run_experiment(experiment, surfaces, ns, objspecs, lb, outpath; nlayouts, verbose)

Sweep one experiment to CSV. Returns the number of cells computed this call.
"""
function run_experiment(experiment::String, surfaces, ns, objspecs, lb, outpath::AbstractString;
                        nlayouts = 30, verbose = true)
    done = completed_keys(outpath)
    ncells = 0
    total = length(surfaces) * length(ns) * length(objspecs) * nlayouts
    t0 = time()
    for (si, s) in enumerate(surfaces)
        buf = NamedTuple[]
        for n in ns, (obj, lam) in objspecs
            lamkey = isnan(lam) ? NaN : lam
            haskey(lb.bank, (n, obj, lamkey)) || continue   # control count skipped as infeasible
            layouts = lb.bank[(n, obj, lamkey)]
            for (li, pa) in enumerate(layouts)
                (s.id, n, obj, lamkey, li) in done && continue
                mask = s.class === :negatives ? pa.negatives : pa.positives
                sum(mask) == n || continue
                geom = layout_geometry(mask, size(s.values)...)
                scores = layout_scores(pa, lb.bounds[n]; lambda = isnan(lam) ? 0.5 : lam)
                for r in sweep_cell(s, mask, geom, scores)
                    push!(buf, merge((surface = s.id, plate = s.meta["plate"],
                                      metric = s.meta["metric"], read = s.meta["read"],
                                      time_h = s.meta["time_h"],
                                      n = n, objective = obj, lambda = lamkey, layout = li), r))
                end
                ncells += 1
            end
        end
        append_rows(outpath, buf)
        if verbose
            el = time() - t0
            @info "[$experiment] surface $si/$(length(surfaces)) ($(s.id)) " *
                  "cells=$ncells elapsed=$(round(el, digits = 1))s"
        end
    end
    return ncells
end
