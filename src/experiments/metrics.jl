# Scoring for predicted error surfaces.
#
# "Truth" throughout is a well's deviation from the plate median, which is what the
# correction actually removes. The absolute OD level is never needed and never known.

"""
    skill(pred, truth) -> Float64

Skill score against the constant baseline, `1 - MSE(pred) / MSE(0)`.

The constant baseline predicts no spatial error at all, i.e. every well gets the
control median. So `0` means the model adds nothing over subtracting that median,
`1` is perfect reconstruction, and negative values mean the model is actively worse.

Being a ratio it is scale-free, so plates whose CV differs by an order of magnitude
can be compared on one axis. Its ceiling is set by the unstructured share of
variance: no smooth interpolator can recover the nugget.
"""
function skill(pred::AbstractVector, truth::AbstractVector)
    denom = mean(abs2, truth)
    denom == 0 && return NaN
    return 1 - mean(abs2, pred .- truth) / denom
end

rmse(pred, truth) = sqrt(mean(abs2, pred .- truth))
mae(pred, truth) = mean(abs, pred .- truth)
bias(pred, truth) = mean(pred .- truth)

"""Spearman rank correlation — agreement of shape, blind to amplitude."""
function spearman(a::AbstractVector, b::AbstractVector)
    length(a) < 3 && return NaN
    ra, rb = tiedrank(a), tiedrank(b)
    s = std(ra) * std(rb)
    return s == 0 ? NaN : cov(ra, rb) / s
end

"""
    coverage(pred, truth, sd; z=1.96) -> Float64

Fraction of wells whose true error falls inside `pred ± z·sd`. Well-calibrated
uncertainty gives ~0.95 at `z=1.96`. This is the metric that silently breaks if a
posterior variance is mistaken for a standard deviation.
"""
function coverage(pred::AbstractVector, truth::AbstractVector, sd::AbstractVector; z = 1.96)
    return mean(abs.(truth .- pred) .<= z .* sd)
end

"""
    variance_decomposition(G) -> NamedTuple

Split a measured surface into row, column, additive row+column, and residual
variance shares, plus a lag-1 nugget estimate and Moran's I.

The nugget share bounds what any smooth interpolator can achieve: it is the part of
well-to-well variance with no spatial structure at the one-well scale, and therefore
sets the plateau height of the accuracy-versus-controls curve.
"""
function variance_decomposition(G::AbstractMatrix)
    A = Float64.(G)
    tot = var(A)
    rowfit = mean(A, dims = 2); colfit = mean(A, dims = 1)
    add = mean(A) .+ (rowfit .- mean(A)) .+ (colfit .- mean(A))
    nr, nc = size(A)

    # lag-1 semivariogram: half the mean squared difference between rook neighbours
    diffs = Float64[]
    for i in 1:nr, j in 1:nc
        i < nr && push!(diffs, (A[i, j] - A[i+1, j])^2)
        j < nc && push!(diffs, (A[i, j] - A[i, j+1])^2)
    end
    nugget = mean(diffs) / 2

    z = vec(A .- mean(A)); n = length(z); W = 0.0; num = 0.0
    for i in 1:nr, j in 1:nc
        k = (j - 1) * nr + i
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            a, b = i + di, j + dj
            (1 <= a <= nr && 1 <= b <= nc) || continue
            W += 1
            num += z[k] * z[(b - 1) * nr + a]
        end
    end

    return (total = tot,
            row_share = 1 - var(A .- rowfit) / tot,
            col_share = 1 - var(A .- colfit) / tot,
            rowcol_share = 1 - var(A .- add) / tot,
            nugget = nugget,
            nugget_share = nugget / tot,
            structured_share = 1 - nugget / tot,
            morans_I = (n / W) * num / sum(abs2, z),
            morans_I_null = -1 / (n - 1))
end
