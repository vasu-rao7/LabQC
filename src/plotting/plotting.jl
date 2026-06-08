

function plot_control_data(data::Matrix{T},plt::PlateArray;kwargs...) where T<: Real 

    size_data_plate_check(data,plt)

    label_types = fill("",size(data)...)

    label_types[findall(x->x==true,runs(plt))] .= "run"
    label_types[findall(x->x==true,plt.negatives)] .= "negative control"
    label_types[findall(x->x==true,plt.positives)] .= "positive control"

    active_labels = label_types[findall(x->x==true,plt.wells)]
    active_data = data[findall(x->x==true,plt.wells)]

    active_labels = CategoricalArray(active_labels)
    levels!(active_labels,["negative control","run","positive control"])

    p = Plots.scatter(active_labels,active_data,color=levelcode.(active_labels),palette=[:red,:black,:blue],legend=false,grid=false)
    xlims!(0.5,3.5,kwargs...)
    return p 
end 



function plot_plate_data(data::Matrix{T},plate::PlateArray;kwargs...) where T<: Real 

    p = plot(plate;kwargs...) 

    heatmap!(data;kwargs...) 

    return p 
end 







