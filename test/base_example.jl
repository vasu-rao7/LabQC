using LabQC , PlateArrays, Distributions, Random 

R,C= 16,24

wells = trues(R,C)

pltval(r,c) = pdf(MvNormal([mean([1,R]),mean([1,C])],[3*R 0 ; 0 4*C]),[r,c])

data = [pltval(i,j) for i in 1:R, j in 1:C]
data = data./maximum(data)

plt  = place_controls(wells,16,16)

p = plot_control_data(data,plt)

p0= plot(plt)

p2 = plot_plate_data(data,plt)

corrected = data_correction(data,plt,"gp_correction")




p3 = plot_plate_data(corrected,plt) 

control_summary(data,plt)

plot(p0,p,p2,p3)


correct, out = LabQC.plateerrorGP(data,plt.negatives,LabQC.absolute;verbose=true)

coords = [transpose(Matrix{Float64}([j i])) for i in 1:R, j in 1:C]

ys = predict_y.((out["gp"],),coords)

yvals = map(x->x[1][1],ys)

avg_neg = mean(data[findall(x->x==true,plt.negatives)])

gp_pred = yvals .+ avg_neg

p4 = plot_plate_data(gp_pred,plt)
