
function datarange(data::Vector{T}) where T<: Real 
    return maximum(data) - minimum(data) 
end 



function control_summary(data::Matrix{T},plate::PlateArray) where T<: Real 

    run = findall(x-> x==true,runs(plate))
    positive = findall(x->x==true , plate.positives) 
    negative = findall(x-> x==true ,plate.negatives)

    idxs = [run,positive,negative]

    data_sets = map(x-> data[x],idxs) 

    df = DataFrame( type = ["runs","positive controls","negative controls"], median = median.(data_sets), mean = mean.(data_sets),sd =std.(data_sets),maximum = maximum.(data_sets),minimum=minimum.(data_sets),range=datarange.(data_sets))

    return df 



end 