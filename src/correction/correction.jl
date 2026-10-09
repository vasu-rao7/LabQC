
"""
    gp_correction(data, plate; verbose=false, passes=:both, kwargs...)

Two-pass GP correction. The first pass fits the negative controls and removes an
additive ("technical") error surface; the second fits the positive controls of the
already-corrected data and applies a multiplicative ("biological") scaling.

`passes` selects which run: `:both` (default), `:technical`, or `:biological`.
An all-control plate has only one class present and must name the pass it wants.

With `verbose=true` the diagnostics of each pass are returned under the keys
`"technical"` and `"biological"`.
"""
function gp_correction(data::Matrix{T}, plate::PlateArray;
                       verbose = false, passes::Symbol = :both, kwargs...) where T <: Real

    passes in (:both, :technical, :biological) ||
        throw(ArgumentError("passes must be :both, :technical or :biological, got :$passes"))

    corrected = convert(Matrix{Float64}, data)
    diagnostics = Dict{String,Any}()

    if passes in (:both, :technical)
        if verbose
            corrected, d = plateerrorGP(corrected, plate.negatives, absolute; verbose = true, kwargs...)
            diagnostics["technical"] = d
        else
            corrected = plateerrorGP(corrected, plate.negatives, absolute; kwargs...)
        end
    end

    if passes in (:both, :biological)
        if verbose
            corrected, d = plateerrorGP(corrected, plate.positives, relative; verbose = true, kwargs...)
            diagnostics["biological"] = d
        else
            corrected = plateerrorGP(corrected, plate.positives, relative; kwargs...)
        end
    end

    return verbose ? (corrected, diagnostics) : corrected
end

function subtract_and_normalize(data::Matrix{T},plate::PlateArray;kwargs...) where T<: Real 

    neg = mean(data[findall(x->x==true,plate.negatives)])
    pos = mean(data[findall(x->x==true,plate.positives)]) 

    return (data .- neg) ./ (pos - neg) 
end 



const data_correction_methods = Dict(
    "gp_correction" => gp_correction,
    "subtract_and_normalize" => subtract_and_normalize
)






"""
    data_correction(data::Matrix{T},plate::PlateArray,method::String;kwargs...) where T <: Real 

Apply a correction method to `data` that have been mapped to a `plate` layout indicating which wells are negative and positive controls 

"""
function data_correction(data::Matrix{T},plate::PlateArray,method::String;kwargs...) where T <: Real 

    size_data_plate_check(data,plate) # check for size errors 

    methodfun  = data_correction_methods[method]

    return methodfun(data,plate;kwargs...)

end
