## MAIN UTILITY FUNCTION

"""
    plateerrorGP(data::Matrix{T}, controls::BitMatrix, correction::Function; kwargs...) where T <: Real

Fit a Gaussian process to the control wells of `data` and use it to estimate, and
remove, position-dependent error at every well.

# Arguments
- `data`: Matrix of raw numeric data
- `controls`: BitMatrix of control locations
- `correction`: method used for correction ([`absolute`](@ref) or [`relative`](@ref))

# Keyword Arguments
- `verbose`: also return a Dict of diagnostics, with keys `"gp"`, `"mean_error"`,
  `"var_f"` (latent posterior variance), `"var_y"` (predictive variance, including
  observation noise) and `"sd_error"` (`sqrt` of `var_f`).
- `ls`, `lb_ls`, `ub_ls`: lengthscale of the isotropic Matern 5/2 kernel and its
  optimisation bounds, in the same units as `pitch`.
- `pitch`: distance between neighbouring wells. Defaults to `1.0`, i.e. well
  indices. Pass `well_pitch(size(data)...)` to fit in mm instead — and change the
  lengthscale bounds to match, since they are expressed in the same units.
- `optimize`: whether to fit kernel hyperparameters. `false` uses `ls`, `sd_guess`
  and `noise` as given.

Only control wells train the model, but every well is corrected.
"""
function plateerrorGP(data::Matrix{T}, controls::BitMatrix, correction::Function;
                      verbose = false, ls = 3.0, lb_ls = 2.0, ub_ls = 6.0,
                      pitch = 1.0, optimize = true, kwargs...) where T <: Real

    ncontrols = sum(controls)
    ncontrols == 0 && error("plateerrorGP: `controls` selects no wells; " *
                            "a plate with no controls of this class cannot be corrected")
    ncontrols < 3 && error("plateerrorGP: only $ncontrols control well(s); at least 3 are " *
                           "needed to fit a lengthscale, signal and noise")
    ncontrols < 6 && @warn "plateerrorGP: fitting $ncontrols control wells; " *
                           "hyperparameters are weakly identified below ~6 controls"

    df = build_df(data, controls)
    avg_control = median(filter(row -> row.Control == 1, df).Value) #calculate median value of set of controls
    df = calculate_controls_error(avg_control, df, correction) #calculate absolute and relative errors and add to dataframe

    # Construct training inputs and outputs from the control wells only
    train_df = filter(row -> row.Control == true, df)
    train_out = Float64.(train_df.CorrectionValue_Train)
    train_in = convert(Matrix{Float64}, copy(transpose(Matrix(train_df[:, [:Column, :Row]])))) .* pitch

    # Fit and optimize GP model
    gp = makeGP(train_in, train_out; lengthscale = ls, lb_lengthscale = lb_ls,
                ub_lengthscale = ub_ls, optimize = optimize, kwargs...)

    # Predict at every well in one call. Looping per well repeats the same
    # triangular solve once per well, which dominates runtime in large sweeps.
    pred_in = convert(Matrix{Float64}, copy(transpose(Matrix(df[:, [:Column, :Row]])))) .* pitch
    error_mu, var_f = predict_f(gp, pred_in)
    var_y = var_f .+ exp(2 * gp.logNoise.value)

    # `absolute` error is an offset; `relative` error is a scaling of the observed value
    adjustby = correction == absolute ? error_mu : df.Value .* error_mu

    df.Error = error_mu
    df.VarF = var_f
    df.VarY = var_y
    df.ErrorSD = sqrt.(max.(var_f, 0.0))
    df.ValueAdjusted = df.Value .- adjustby

    d_corrected = zeros(Float64, size(data)...)
    d_mean_error = zeros(Float64, size(data)...)
    d_var_f = zeros(Float64, size(data)...)
    d_var_y = zeros(Float64, size(data)...)
    d_sd_error = zeros(Float64, size(data)...)

    for row in eachrow(df)
        d_corrected[row.Row, row.Column] = row.ValueAdjusted
        d_mean_error[row.Row, row.Column] = row.Error
        d_var_f[row.Row, row.Column] = row.VarF
        d_var_y[row.Row, row.Column] = row.VarY
        d_sd_error[row.Row, row.Column] = row.ErrorSD
    end

    if verbose == true
        return d_corrected, Dict("gp" => gp,
                                 "mean_error" => d_mean_error,
                                 "var_f" => d_var_f,
                                 "var_y" => d_var_y,
                                 "sd_error" => d_sd_error)
    end
    return d_corrected

end

"""
    makeGP(train_in, train_out; kwargs...)

Zero-mean GP with an isotropic Matern 5/2 kernel over well positions.

`optimize = false` skips hyperparameter fitting and uses `lengthscale`, `sd_guess`
and `noise` as supplied. No catch-all keyword argument is accepted: a misspelled
hyperparameter must raise a `MethodError` rather than being silently discarded,
because `GaussianProcesses.optimize!` drops unrecognised keywords without warning
when kernel bounds are given.
"""
function makeGP(train_in, train_out; lengthscale = 3.0, sd_guess = signal_guess(train_out),
                lb_lengthscale = 2.0, ub_lengthscale = 6.0,
                lb_std = sd_guess / 5, ub_std = sd_guess * 3,
                noise = -2.0, optimize = true)
    mZero = MeanZero()

    kern = Matern(5/2, log(lengthscale), log(sd_guess))

    # Fit a GP model to training data
    gp = GP(train_in, train_out, mZero, kern, noise)

    optimize || return gp

    # Define kernel bounds for GP hyperparam optimization
    kb = [
            [log(lb_lengthscale), log(lb_std)],   # lower bounds
            [log(ub_lengthscale), log(ub_std)]    # upper bound
        ]

    # Optimize GP
    GaussianProcesses.optimize!(gp; domean = false, kernbounds = kb)

    return gp
end

"""
    signal_guess(train_out) -> Float64

Starting signal standard deviation. `std` is degenerate for the small control sets
this model is fit to — undefined for one well, and exactly zero when every control
reads the same — and `log(0)` produces a kernel full of `NaN`, so fall back to a
small positive value in those cases.
"""
function signal_guess(train_out)
    length(train_out) < 2 && return 1.0
    s = std(train_out)
    return (!isfinite(s) || s <= 0) ? 1.0 : s
end

function build_df(data::Matrix, controls::BitMatrix)
    size(data) == size(controls) ? nothing : error("size of data must be the size of controls")


    # Initialize vectors to store the values, row indices, and column indices
    values = Float64[]
    row_indices = Int[]
    col_indices = Int[]


    # Iterate over the matrix and fill the vectors
    for idx in CartesianIndices(data)
            push!(values, data[idx])
            push!(row_indices, idx[1])
            push!(col_indices, idx[2])

    end

    # Create a DataFrame from the vectors
    df = DataFrame(Value = values, Row = row_indices, Column = col_indices, Control = vec(controls))
    return df
end

## Calculate the error (either absolute or relative and add these values to dataframe)
function calculate_controls_error(avg_control, df::DataFrame,correction::Function)
    avg_val =Float64[]
    corr_val = Float64[]

    for i in eachrow(df)
        if i.Control==1
            push!(avg_val, avg_control)
        else
            push!(avg_val, NaN64)
        end
        ybar = correction(i.Value, avg_val[end])
        push!(corr_val, ybar)

    end
    df.AverageValue = avg_val
    df.CorrectionValue_Train = corr_val


    return df
end


absolute(x,y) = x-y 
relative(x,y) = (x-y)/x 
