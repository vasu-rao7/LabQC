# Ingest tests. The fixtures are the project's own raw exports, chosen because each
# one exercises a hazard that silently corrupts data if the reader gets it wrong.

using Test, LabQC, Statistics, Dates
using LabQC: tonum, totime_h, unwrap_hours!, grid_hash

const RAW = joinpath(@__DIR__, "..", "raw_data")

@testset "cell parsing" begin
    @test tonum(0.437) === 0.437
    @test tonum("0.437") === 0.437
    @test tonum(missing) === missing
    # BioTek writes these for saturated / failed reads; they must not become numbers
    @test tonum("?????") === missing
    @test tonum("OVRFLW") === missing
    @test tonum("") === missing
    @test tonum("not a number") === missing
end

@testset "elapsed-time parsing" begin
    @test totime_h(Dates.Time(0, 0, 12)) ≈ 12 / 3600
    @test totime_h(Dates.Time(13, 10, 12)) ≈ 13 + 10 / 60 + 12 / 3600
    # a read past 24 h is typed as a DateTime on the Excel epoch, not a Time
    @test totime_h(Dates.DateTime("1899-12-31T00:00:10")) ≈ 24 + 10 / 3600
    @test totime_h(1.0) ≈ 24.0
    @test totime_h(missing) === missing
    # Dates.Time wraps at midnight; elapsed time must stay monotonic
    @test unwrap_hours!([23.5, 23.75, 0.003]) ≈ [23.5, 23.75, 24.003]
end

@testset "kinetic: header must not be trusted over the data" begin
    k = read_kinetic(joinpath(RAW, "092424_neg_control_plate1.xlsx"))
    # the file declares "Interval 0:00:15, 97 Reads" (24 h) and is wrong on both counts
    @test occursin("97 Reads", k.meta["Discontinuous Kinetics"])
    @test k.nreads == 80                     # not 97
    @test k.interval_min ≈ 10.0 atol = 0.5   # not 15
    @test k.blank_rows == 17                 # padded to the declared count
    @test k.times_h[end] ≈ 13.17 atol = 0.01
    @test size(k.values) == (8, 12, 80)
    @test k.n_missing == 0
    @test issorted(k.times_h)
end

@testset "kinetic: a full-length run reads all of its data" begin
    k = read_kinetic(joinpath(RAW, "2026_08_21_VR_PosControl.xlsx"))
    @test k.nreads == 97
    @test k.interval_min ≈ 15.0 atol = 0.5
    @test k.blank_rows == 0
    @test k.times_h[end] ≈ 24.0 atol = 0.01   # the 24 h read is not silently dropped
    @test issorted(k.times_h)
    # peak is the plate mean maximum, ~8.5 h in
    m = [mean(skipmissing(k.values[:, :, t])) for t in 1:k.nreads]
    @test k.times_h[argmax(m)] ≈ 8.5 atol = 0.3
end

@testset "transient is flagged, not dropped" begin
    k = read_kinetic(joinpath(RAW, "092424_neg_control_plate1.xlsx"))
    tr = transient_reads(k)
    @test first(tr) == 1
    @test 2 <= last(tr) <= 10          # settling occupies the first few reads only
    @test size(k.values, 3) == k.nreads # nothing was removed
end

@testset "endpoint grids and duplicate detection" begin
    files = filter(f -> endswith(f, ".xlsx"), readdir(joinpath(RAW, "2026-06-04")))
    grids = reduce(vcat, [read_endpoint_grids(joinpath(RAW, "2026-06-04", f)) for f in files])
    @test length(grids) == 14                       # sheet-level grids
    @test length(unique(g.hash for g in grids)) == 8  # genuinely distinct surfaces
    @test all(g -> size(g.values) == (8, 12), grids)

    by_src = Dict(string(g.source, "|", g.sheet) => g.hash for g in grids)
    # known duplicates: exporting the same plate twice must not look like replication
    @test by_src["plate1_read2.xlsx|Sheet1"] ==
          by_src["plate3_read2.xlsx|Plate 2 - Sheet1"]
    @test by_src["plate1_read3.xlsx|Plate 2 - Sheet1"] ==
          by_src["plate3_read3.xlsx|Plate 2 - Sheet1"]
    # and genuinely different reads must not collide
    @test by_src["plate1_read1.xlsx|Plate 1 - Sheet1"] !=
          by_src["plate1_read2.xlsx|Sheet1"]

    od = filter(g -> get(g.meta, "Read", "") == "OD600", grids)
    @test length(unique(g.hash for g in od)) == 5
end

@testset "the gas gradient is present and column-structured" begin
    g = only(filter(x -> x.hash == first([y.hash for y in
            read_endpoint_grids(joinpath(RAW, "2026-06-04", "plate1_read1.xlsx"))]),
        read_endpoint_grids(joinpath(RAW, "2026-06-04", "plate1_read1.xlsx"))))
    G = Float64.(g.values)
    colmeans = vec(mean(G, dims = 1))
    rowmeans = vec(mean(G, dims = 2))
    # column 12 sits nearest the fan and is strongly suppressed relative to the plateau
    @test colmeans[12] < 0.6 * mean(colmeans[3:9])
    # column position explains far more variance than row position
    varexpl(M, d) = 1 - var(M .- mean(M, dims = d)) / var(M)
    @test varexpl(G, 1) > 0.5      # removing column means kills most of the variance
    @test varexpl(G, 2) < 0.3
end
