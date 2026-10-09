# Readers for BioTek plate-reader .xlsx exports.
#
# Two export shapes appear in this project's raw data:
#   * endpoint   - an 8x12 grid, row letters down column c, values in the 12 columns to the right
#   * kinetic    - a wide table, one row per read, wells A1..H12 spread across the columns
#
# Neither shape is at a fixed cell address, so both readers anchor on a token
# ("A" for endpoint, "A1" for kinetic) and validate that exactly 96 wells were found.
#
# Three hazards these readers exist to handle, all observed in raw_data/:
#   1. The declared read count can be wrong. 092424_neg_control_plate1 announces
#      "97 Reads, Interval 0:00:15" but holds 80 reads at 10 min, padded with 17
#      blank rows whose Time cell is 0. Never trust the header; stop at the padding.
#   2. Saturated or failed reads are written as "?????" or "OVRFLW", not as numbers.
#   3. A single file may carry several sheets, and different files may carry
#      byte-identical grids. Hash every grid so duplicates cannot be counted twice.

const WELL_RE = r"^([A-H])(1[0-2]|[1-9])$"
const SENTINELS = ("?????", "OVRFLW", "OVER", "")

"""Parse a spreadsheet cell into an OD value, mapping BioTek sentinels to `missing`."""
function tonum(x)
    x === missing && return missing
    x isa Real && return Float64(x)
    s = strip(string(x))
    (s in SENTINELS) && return missing
    v = tryparse(Float64, s)
    return v === nothing ? missing : v
end

"""
Elapsed time of a read, in hours.

XLSX.jl types BioTek's elapsed-time column as `Dates.Time`, which wraps at 24 h;
`unwrap_hours!` repairs that. A plain number is an Excel day fraction.
"""
const EXCEL_EPOCH = Dates.DateTime(1899, 12, 30)

function totime_h(x)
    x === missing && return missing
    x isa Dates.Time && return Dates.value(x) / 3.6e12       # ns -> h
    # a read past the 24 h mark is typed as a DateTime on the Excel epoch
    x isa Dates.DateTime && return Dates.value(x - EXCEL_EPOCH) / 3.6e6   # ms -> h
    x isa Dates.Date && return Dates.value(x - Dates.Date(EXCEL_EPOCH)) * 24.0
    x isa Real && return Float64(x) * 24.0                   # day fraction -> h
    v = tryparse(Float64, strip(string(x)))
    return v === nothing ? missing : v * 24.0
end

"""Restore monotonicity to elapsed times that wrapped at the 24 h boundary."""
function unwrap_hours!(t::Vector{Float64})
    off = 0.0
    for k in 2:length(t)
        t[k] + off < t[k-1] && (off += 24.0)
        t[k] += off
    end
    return t
end

"""Stable content hash of a plate grid, used to detect duplicate exports."""
function grid_hash(G::AbstractMatrix)
    io = IOBuffer()
    for v in G
        print(io, v === missing ? "NA" : string(round(Float64(v), digits = 6)), ",")
    end
    return bytes2hex(SHA.sha256(take!(io)))[1:12]
end

"""Pull the `Key => value` metadata rows BioTek writes above the data block."""
function sheet_meta(D::AbstractMatrix)
    keys_wanted = ("Read", "Plate Number", "Reader Type:", "Plate Type", "Date",
                   "Software Version", "Experiment File Path:", "Discontinuous Kinetics")
    meta = Dict{String,Any}()
    for r in 1:size(D, 1), c in 1:min(size(D, 2), 3)
        k = D[r, c]
        k isa AbstractString || continue
        if strip(k) in keys_wanted && c + 1 <= size(D, 2)
            v = D[r, c+1]
            v === missing || (meta[strip(k)] = v)
        end
    end
    return meta
end

sheets_of(path) = [(n, XLSX.getdata(XLSX.readxlsx(path)[n])) for n in XLSX.sheetnames(XLSX.readxlsx(path))]

"""
    read_endpoint_grids(path) -> Vector{NamedTuple}

Every 8x12 endpoint grid in the workbook, one entry per grid, each carrying its
source sheet, a content hash, and the sheet metadata.
"""
function read_endpoint_grids(path::AbstractString)
    out = NamedTuple[]
    for (sheet, D) in sheets_of(path)
        meta = sheet_meta(D)
        nr, nc = size(D)
        for r in 1:nr, c in 1:nc
            # anchor: a cell reading "A" with "B".."H" directly beneath it
            D[r, c] isa AbstractString && strip(string(D[r, c])) == "A" || continue
            r + 7 <= nr && c + 12 <= nc || continue
            all(i -> D[r+i-1, c] isa AbstractString &&
                     strip(string(D[r+i-1, c])) == string("ABCDEFGH"[i]), 1:8) || continue
            G = Matrix{Union{Float64,Missing}}(undef, 8, 12)
            for i in 1:8, j in 1:12
                G[i, j] = tonum(D[r+i-1, c+j])
            end
            count(!ismissing, G) >= 90 || continue
            push!(out, (values = G, source = basename(path), sheet = sheet,
                        hash = grid_hash(G), meta = meta))
        end
    end
    return out
end

"""
    read_kinetic(path) -> NamedTuple

A kinetic export as `values` (8x12xT), `times_h`, and metadata. Stops at the first
all-blank read row so header-declared read counts are never trusted, and derives the
read interval from the timestamps rather than from the declared protocol.
"""
function read_kinetic(path::AbstractString)
    for (sheet, D) in sheets_of(path)
        nr, nc = size(D)
        hdr = findfirst(r -> any(c -> D[r, c] isa AbstractString &&
                                 strip(string(D[r, c])) == "A1", 1:nc), 1:nr)
        hdr === nothing && continue

        wellcol = Dict{String,Int}()
        timecol = nothing
        for c in 1:nc
            v = D[hdr, c]
            v isa AbstractString || continue
            s = strip(string(v))
            occursin(WELL_RE, s) && (wellcol[s] = c)
            s == "Time" && (timecol = c)
        end
        length(wellcol) == 96 || error("$(basename(path)): found $(length(wellcol)) wells, expected 96")
        timecol === nothing && error("$(basename(path)): no Time column in header row $hdr")

        times = Float64[]
        frames = Matrix{Union{Float64,Missing}}[]
        stopped_at = nr + 1
        for r in (hdr+1):nr
            row = Matrix{Union{Float64,Missing}}(undef, 8, 12)
            nfilled = 0
            for (w, c) in wellcol
                m = match(WELL_RE, w)
                v = tonum(D[r, c])
                row[Int(m[1][1]) - Int('A') + 1, parse(Int, m[2])] = v
                ismissing(v) || (nfilled += 1)
            end
            if nfilled == 0
                isempty(frames) && continue      # blanks before the block starts
                stopped_at = r                   # padding: the data block ends here
                break
            end
            t = totime_h(D[r, timecol])
            ismissing(t) && continue
            push!(times, Float64(t))
            push!(frames, row)
        end

        # count only the run of padded blank rows directly after the data block,
        # so trailing sections (BioTek's "Results" table) are not miscounted
        blanks = 0
        r = stopped_at
        while r <= nr && !ismissing(totime_h(D[r, timecol])) &&
              all(c -> ismissing(tonum(D[r, c])), values(wellcol))
            blanks += 1
            r += 1
        end

        isempty(frames) && error("$(basename(path)): no data rows found")
        unwrap_hours!(times)

        V = Array{Union{Float64,Missing}}(undef, 8, 12, length(frames))
        for (k, f) in enumerate(frames); V[:, :, k] = f; end
        ivals = length(times) > 1 ? diff(times) .* 60 : Float64[]

        return (values = V, times_h = times, source = basename(path), sheet = sheet,
                nreads = length(times),
                interval_min = isempty(ivals) ? NaN : median(ivals),
                blank_rows = blanks,
                n_missing = count(ismissing, V),
                meta = sheet_meta(D))
    end
    error("$(basename(path)): no kinetic block found")
end

"""
    transient_reads(k; frac = 0.9) -> UnitRange

Reads belonging to the start-up transient, where the plate mean is still falling
from its initial settling value. Returned as metadata; nothing is dropped.
"""
function transient_reads(k; frac = 0.9)
    m = [mean(skipmissing(k.values[:, :, t])) for t in 1:k.nreads]
    lo = argmin(m)
    last_t = 0
    for t in 1:lo
        m[t] > m[lo] / frac && (last_t = t)
    end
    return 1:last_t
end
