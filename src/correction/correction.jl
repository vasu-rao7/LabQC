
function gp_correction(data::Matrix{T},plate::PlateArray;kwargs...) where T<: Real 

    neg_corrected = plateerrorGP(data,plate.negatives, absolute;kwargs...)
    pos_corrected = plateerrorGP(neg_corrected,plate.positives,relative;kwargs...)

    return pos_corrected

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