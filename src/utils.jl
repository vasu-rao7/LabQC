function size_data_plate_check(data::Matrix{T},plate::PlateArray) where T<: Real 
    size(data) == size(plate) || error("Size of data matrix ($(size(data)) must equal size of PlateArray ($(size(plate))))")
    return nothing 
end 