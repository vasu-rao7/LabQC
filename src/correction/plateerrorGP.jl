## MAIN UTILITY FUNCTION

"""
    plateerrorGP(data::Matrix{T}, controls::BitMatrix, correction::Function) where T <: Real

# Arguments 
- data: Matrix of raw numeric data
- controls: BitMatrix of control locations 
- correction: method used for correction (absolute or relative)
- verbose: returns a dictionary containing the gp object, a matrix of the GP's mean error predictions, and a matrix of the GP's standard deviations. The dictionary has keys "gp", "mean_error", and "sd_error"
- ls: lengthscale of isotropic Mater52 GP kernel
- lb_ls: lower bound of lengthscale
- ub_ls: upper bound of lengthscale

"""
function plateerrorGP(data::Matrix{T}, controls::BitMatrix, correction::Function; verbose=false, ls=3.0, lb_ls=2.0, ub_ls=6.0,kwargs...) where T <: Real

    df = build_df(data, controls)
    # df = DataFrame(Value = values, Row = row_indices, Column = col_indices, Control = controls)
    avg_control = median(filter(row -> row.Control==1, df).Value) #calculate median value of set of controls
    df = calculate_controls_error(avg_control, df, correction) #calculate absolute and relative errors and add to dataframe
    # Construct Training inputs and outputs
    train_out = Float64[]
    for row in eachrow(df)
        if row.Control
            push!(train_out, row.CorrectionValue_Train)
        end
    end
    train_filter_df = filter(row -> row.Control == true, df)
    train_in = copy(transpose(Matrix(train_filter_df[:, [:Column, :Row]])))
    train_in = convert(Matrix{Float64}, train_in)

    # Fit and optimize GP model
    gp = makeGP(train_in, train_out; lengthscale=ls, lb_lengthscale=lb_ls, ub_lengthscale=ub_ls,kwargs...)

    # Get GP predictions for all wells and add these values to dataframe
    error_μ = Float64[]
    error_sd = Float64[]
    val_adjusted = Float64[]

    for i in eachrow(df)
        h = hcat(i.Column, i.Row)
        h = convert(Matrix{Float64}, h)
        pred_pt=copy(transpose(h))
        pred_error_μ, pred_error_sd = predict_y(gp, pred_pt)
        if correction == absolute
            adjustby = pred_error_μ[1]
        else
            adjustby = i.Value * pred_error_μ[1]
        end
        val_adj = i.Value - adjustby
        push!(error_μ, pred_error_μ[1])
        push!(error_sd, pred_error_sd[1])
        push!(val_adjusted, val_adj)
    end

    df.Error = error_μ
    df.ErrorSD = error_sd
    df.ValueAdjusted = val_adjusted

    extra_outputs=Dict()

    d_corrected = zeros(Float64, size(data)[1], size(data)[2])
    d_mean_error = zeros(Float64, size(data)[1], size(data)[2])
    d_sd_error = zeros(Float64, size(data)[1], size(data)[2])
    
    for row in eachrow(df)
        d_corrected[row.Row, row.Column] = row.ValueAdjusted
        d_mean_error[row.Row, row.Column] = row.Error
        d_sd_error[row.Row, row.Column] = row.ErrorSD
    end
    if verbose == true
        extra_outputs["gp"]=gp
        extra_outputs["mean_error"] = d_mean_error
        extra_outputs["sd_error"] = d_sd_error
        return d_corrected, extra_outputs
    end
    return d_corrected
    
end

function makeGP(train_in, train_out; lengthscale=3.0, sd_guess=std(train_out), lb_lengthscale=2.0, ub_lengthscale=6.0, lb_std = std(train_out)/5, ub_std=std(train_out)*3,noise=-2.0,kwargs...)
    mZero = MeanZero()


    kern = Matern(5/2, log(lengthscale), log(sd_guess))

    # Fit a GP model to training data
    gp = GP(train_in,train_out,mZero,kern,noise)

    # Define kernel bounds for GP hyperparam optimization
    kb = [
            [log(lb_lengthscale), log(lb_std)],   # lower bounds
            [log(ub_lengthscale), log(ub_std)]    # upper bound
        ]

    # Optimize GP
    GaussianProcesses.optimize!(gp; domean=false, kernbounds = kb,kwargs...)

    return gp
end

function build_df(data::Matrix, controls::BitMatrix)
    size(data) == size(controls) ? nothing : error("size of data must be the size of controls")


    # Initialize vectors to store the values, row indices, and column indices
    values = Float64[]
    row_indices = Int[]
    col_indices = Int[]
    is_control = Bool[]


    # Iterate over the matrix and fill the vectors
    for idx in CartesianIndices(data)
            push!(values, data[idx])
            push!(is_control,controls[idx])
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






