using LabQC, PlateArrays ,Distributions 

Random.seed!(123345)

R,C = 8,12 

wells = trues(R,C) 

pos_dist = Normal(1.5,0.15)

neg_dist = Normal(0.12,0.06)
rundist = Uniform(0.25,1.6)

data = zeros(R,C)
n_pos,n_neg = 8,8 

nruns = R*C - n_pos - n_neg 
plt  = place_controls(wells,n_pos,n_neg)

posvals = rand(pos_dist,n_pos)
negvals = rand(neg_dist,n_neg)
runvals = rand(rundist,nruns)

data[findall(x->x==true,runs(plt))] .= runvals 
data[findall(x->x==true,plt.positives)] .= posvals
data[findall(x->x==true,plt.negatives)] .= negvals

p0 = plot(plt) 

p1 = plot_plate_data(data,plt)

p2 = plot_control_data(data,plt) 
hline!([0.25],color="red")


savefig(plot(p0,p1,p2,layout=grid(1,3, widths=(3.5/10,4.5/10,2/10)),size=(1600,400)),"/Users/BDavid/Desktop/auto_sci_controls.svg")


using LabQC, PlateArrays ,Distributions 

Random.seed!(123345)

R,C = 8,12 

wells = trues(R,C) 

pos_dist = Normal(1.5,0.15)

neg_dist = Normal(0.12,0.06)
rundist = Uniform(0.25,1.6)

data = zeros(R,C)
n_pos,n_neg = 8,8 

nruns = R*C - n_pos - n_neg 
plt  = place_controls(wells,n_pos,n_neg)

posvals = rand(pos_dist,n_pos)
negvals = rand(neg_dist,n_neg)
 
negvals[1:2] .= 0.3 .+ 0.3* rand(2)
runvals = rand(rundist,nruns)

data[findall(x->x==true,runs(plt))] .= runvals 
data[findall(x->x==true,plt.positives)] .= posvals
data[findall(x->x==true,plt.negatives)] .= negvals

p0 = plot(plt) 

p1 = plot_plate_data(data,plt)

p2 = plot_control_data(data,plt) 
hline!([0.25],color="red")


savefig(plot(p0,p1,p2,layout=grid(1,3, widths=(3.5/10,4.5/10,2/10)),size=(1600,400)),"/Users/BDavid/Desktop/auto_sci_controls_bad.svg")
