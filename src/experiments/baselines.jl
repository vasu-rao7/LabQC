# Alternative error-surface predictors, for comparison against the GP.
#
# Each takes the measured plate, a control mask, and the correction whose scale the
# errors live on, and returns a predicted error at every well. All are fit ONLY on
# the control wells, so they face exactly the problem the GP faces, on the same scale.

"""
    truth_surface(G, mask, correction) -> Matrix

The error each model is asked to predict, on the same scale the correction removes:
deviation from the control median for [`absolute`](@ref), and the fractional version
of that for [`relative`](@ref). This is also the GP's training target at the controls.
"""
function truth_surface(G::AbstractMatrix, mask::BitMatrix, correction::Function)
    A = Float64.(G)
    m = median(A[mask])
    return correction === absolute ? A .- m : (A .- m) ./ A
end

"""
    error_target(G, mask, correction) -> Matrix

The quantity every model predicts, defined at all wells but observable to a model
fit on `mask` only at the control entries. There is deliberately no default
`correction`: predicting on the wrong scale is the bug this signature exists to stop.
"""
error_target(G::AbstractMatrix, mask::BitMatrix, correction::Function) =
    truth_surface(G, mask, correction)

"""Predict no spatial error anywhere. The null model, and the skill-score denominator."""
predict_constant(G, mask, correction) = zeros(Float64, size(G))

"""
    predict_rowcol(G, mask, correction) -> Matrix

Additive row + column model fit to the control wells by minimum-norm least squares.

Note this model is *underdetermined* at small control counts: an 8x12 plate has
19 free effects (8 rows + 12 columns, less one for the shared intercept) against
only `n` observations. The minimum-norm solution is the natural choice there, but
the model cannot be expected to shine until `n` approaches 20. Reported anyway,
because additive row/column correction is the incumbent in the plate-QC literature.
"""
function predict_rowcol(G::AbstractMatrix, mask::BitMatrix, correction::Function)
    nr, nc = size(G)
    idx = findall(mask)
    t = error_target(G, mask, correction)
    y = [t[i] for i in idx]
    X = zeros(Float64, length(idx), nr + nc)
    for (k, ci) in enumerate(idx)
        X[k, ci[1]] = 1.0
        X[k, nr + ci[2]] = 1.0
    end
    beta = pinv(X) * y                      # minimum-norm solution when rank-deficient
    return [beta[i] + beta[nr + j] for i in 1:nr, j in 1:nc]
end

"""Predict each well's error as that of its nearest control well."""
function predict_nearest(G::AbstractMatrix, mask::BitMatrix, correction::Function)
    nr, nc = size(G)
    t = error_target(G, mask, correction)
    idx = findall(mask)
    out = zeros(Float64, nr, nc)
    for i in 1:nr, j in 1:nc
        best, bd = 0.0, Inf
        for ci in idx
            d = (ci[1] - i)^2 + (ci[2] - j)^2
            d < bd && (bd = d; best = t[ci])
        end
        out[i, j] = best
    end
    return out
end

"""Inverse-distance-weighted interpolation of control-well error. Smooth, but not probabilistic."""
function predict_idw(G::AbstractMatrix, mask::BitMatrix, correction::Function; power = 2.0)
    nr, nc = size(G)
    t = error_target(G, mask, correction)
    idx = findall(mask)
    out = zeros(Float64, nr, nc)
    for i in 1:nr, j in 1:nc
        num, den, exact = 0.0, 0.0, NaN
        for ci in idx
            d2 = (ci[1] - i)^2 + (ci[2] - j)^2
            if d2 == 0
                exact = t[ci]; break
            end
            w = d2^(-power / 2)
            num += w * t[ci]; den += w
        end
        out[i, j] = isnan(exact) ? num / den : exact
    end
    return out
end

const baseline_methods = Dict(
    "constant" => predict_constant,
    "rowcol"   => predict_rowcol,
    "nearest"  => predict_nearest,
    "idw"      => predict_idw,
)
