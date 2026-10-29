# GP correction tests. These pin the fixes made for the experiment harness:
# a correct posterior sd, guards for degenerate control sets, unit-independent
# geometry, and single-pass correction for all-control plates.

using Test, LabQC, PlateArrays, Statistics, Random
using LabQC: plateerrorGP, absolute, relative, signal_guess

const R, C = 8, 12
bump(r, c) = exp(-((r - 4.5)^2 / 12 + (c - 6.5)^2 / 28))
const DATA = [0.2 + 0.5 * bump(i, j) for i in 1:R, j in 1:C]

function split_masks(seed = 11)
    Random.seed!(seed)
    neg = falses(R, C); pos = falses(R, C)
    for (k, idx) in enumerate(shuffle(vec(CartesianIndices(DATA)))[1:24])
        (k <= 12 ? neg : pos)[idx] = true
    end
    return neg, pos
end

@testset "posterior spread is a standard deviation, not a variance" begin
    neg, _ = split_masks()
    _, d = plateerrorGP(DATA, neg, absolute; verbose = true)
    # the historical bug: predict_y's second return is a VARIANCE, stored as "sd_error"
    @test d["sd_error"] ≈ sqrt.(d["var_f"])
    @test all(d["var_y"] .>= d["var_f"])        # predictive variance includes noise
    @test all(d["var_f"] .>= 0)
    # a variance mistaken for an sd makes intervals far too narrow
    @test mean(d["sd_error"]) > 3 * mean(d["var_f"])
end

@testset "degenerate control sets fail loudly" begin
    neg, _ = split_masks()
    @test_throws ErrorException plateerrorGP(DATA, falses(R, C), absolute)
    two = falses(R, C); two[1, 1] = true; two[8, 12] = true
    @test_throws ErrorException plateerrorGP(DATA, two, absolute)
    # a misspelled hyperparameter must not be silently discarded
    @test_throws MethodError plateerrorGP(DATA, neg, absolute; lenghtscale = 27.0)
end

@testset "signal_guess survives degenerate inputs" begin
    @test signal_guess([1.0]) > 0                 # std undefined for one point
    @test signal_guess(fill(0.5, 8)) > 0          # std exactly zero
    @test signal_guess([1.0, 2.0, 3.0]) ≈ std([1.0, 2.0, 3.0])
    @test isfinite(log(signal_guess(fill(0.5, 8))))
end

@testset "geometry: mm and well-index fits agree when bounds scale together" begin
    neg, _ = split_masks()
    @test well_pitch(8, 12) == 9.0
    @test well_pitch(16, 24) == 4.5
    @test_throws ErrorException well_pitch(7, 13)
    a = plateerrorGP(DATA, neg, absolute; ls = 3.0, lb_ls = 2.0, ub_ls = 6.0)
    p = well_pitch(R, C)
    b = plateerrorGP(DATA, neg, absolute; pitch = p, ls = 3p, lb_ls = 2p, ub_ls = 6p)
    @test a ≈ b
end

@testset "all-control plates can run a single pass" begin
    allneg = PlateArray(trues(R, C), falses(R, C), trues(R, C))
    # the two-pass default cannot work when a control class is absent
    @test_throws ErrorException gp_correction(DATA, allneg)
    out = gp_correction(DATA, allneg; passes = :technical)
    @test size(out) == size(DATA)
    @test all(isfinite, out)
    @test_throws ArgumentError gp_correction(DATA, allneg; passes = :nonsense)
end

@testset "verbose returns per-pass diagnostics" begin
    neg, pos = split_masks()
    plate = PlateArray(trues(R, C), pos, neg)
    # this previously threw: pass 1's tuple was fed to pass 2 as the data matrix
    out, diag = gp_correction(DATA, plate; verbose = true)
    @test size(out) == size(DATA)
    @test Set(keys(diag)) == Set(["technical", "biological"])
    @test haskey(diag["technical"], "var_f")
    @test size(diag["biological"]["mean_error"]) == size(DATA)
    out2, diag2 = gp_correction(DATA, plate; verbose = true, passes = :technical)
    @test collect(keys(diag2)) == ["technical"]
    @test out2 ≈ plateerrorGP(DATA, neg, absolute)
end

@testset "correction removes the injected spatial structure" begin
    neg, pos = split_masks()
    plate = PlateArray(trues(R, C), pos, neg)
    corrected = gp_correction(DATA, plate; passes = :technical)
    @test std(corrected) < 0.5 * std(DATA)
end

@testset "default kernel is unchanged by the kernel option" begin
    # golden values captured at commit 6755b12, before `kernel` existed
    golden = Dict{Tuple{String,String,Int,Int},Float64}()
    for line in eachline(joinpath(@__DIR__, "golden", "plateerrorGP_iso.txt"))
        startswith(line, "#") && continue
        case, key, i, j, v = split(line, '\t')
        golden[(case, key, parse(Int, i), parse(Int, j))] = parse(Float64, v)
    end
    neg, pos = split_masks()
    for (case, mask, corr) in (("neg_absolute", neg, absolute), ("pos_relative", pos, relative))
        corrected, d = plateerrorGP(DATA, mask, corr; verbose = true)
        for (key, M) in (("corrected", corrected), ("mean_error", d["mean_error"]),
                         ("var_f", d["var_f"]), ("var_y", d["var_y"]))
            @test all(M[i, j] === golden[(case, key, i, j)] for i in 1:R, j in 1:C)
        end
    end
end

@testset "per-axis (ARD) kernel" begin
    neg, _ = split_masks()
    _, d = plateerrorGP(DATA, neg, absolute; verbose = true, kernel = :ard, ls = (4.0, 3.0))
    gp = d["gp"]
    @test length(LabQC.lengthscales(gp)) == 2
    @test all(2.0 - 1e-6 .<= LabQC.lengthscales(gp) .<= 6.0 + 1e-6)
    _, d_iso = plateerrorGP(DATA, neg, absolute; verbose = true)
    @test length(LabQC.lengthscales(d_iso["gp"])) == 1
    @test LabQC.signal_sd(d_iso["gp"]) ≈ sqrt(d_iso["gp"].kernel.σ2)
    @test_throws ArgumentError plateerrorGP(DATA, neg, absolute; kernel = :nonsense)
end
