using JuMP
using Gurobi
using LinearAlgebra
using CSV
using DataFrames

include("GML_struct.jl")
include("GML.jl")

# process data
# reading data
Bank1_raw=CSV.read("data/demand_and_solar_allendale_bank139_aug2016.csv");
Bank2_raw=CSV.read("data/demand_and_solar_allendale_bank239_aug2016.csv");
Bank3_raw=CSV.read("data/demand_and_solar_blooming_grove_bank276_aug2016.csv");
# initial parameters
Bank1_Pd_raw_5th=zeros(288,4)
Bank2_Pd_raw_5th=zeros(288,4)
Bank3_Pd_raw_5th=zeros(288,4)
Bank1_Qd_raw_5th=zeros(288,4)
Bank2_Qd_raw_5th=zeros(288,4)
Bank3_Qd_raw_5th=zeros(288,4)
Bank1_Sg_raw_5th=zeros(288,4)
Bank2_Sg_raw_5th=zeros(288,4)
Bank3_Sg_raw_5th=zeros(288,4)

# get data
# Here we read the data of Feb 1st.
# Reading the data of other day just shift the k in the for loop.
for i=1:288
    k=(i-1)*5+1;
    # bank 1
    Bank1_Pd_raw_5th[i,:]=[Bank1_raw[k,4],Bank1_raw[k,9],Bank1_raw[k,14],
    Bank1_raw[k,19]]
    Bank1_Qd_raw_5th[i,:]=Bank1_raw[k,3]*[Bank1_raw[k,4],Bank1_raw[k,9],
    Bank1_raw[k,14], Bank1_raw[k,19]]
    Bank1_Sg_raw_5th[i,:]=[1080, 1155, 1790, 1775]*Bank1_raw[k,24]
    # bank 2
    Bank2_Pd_raw_5th[i,:]=[Bank2_raw[k,4],Bank2_raw[k,9],Bank2_raw[k,14],
    Bank2_raw[k,19]]
    Bank2_Qd_raw_5th[i,:]=Bank2_raw[k,3]*[Bank2_raw[k,4],Bank2_raw[k,9],
    Bank2_raw[k,14], Bank2_raw[k,19]]
    Bank3_Sg_raw_5th[i,:]=[1080, 1155, 1790, 1775]*Bank2_raw[k,24]
    # bank 3
    Bank3_Pd_raw_5th[i,:]=[Bank3_raw[k,4],Bank3_raw[k,9],Bank3_raw[k,14],
    Bank3_raw[k,19]]
    Bank3_Qd_raw_5th[i,:]=Bank3_raw[k,3]*[Bank3_raw[k,4],Bank3_raw[k,9],
    Bank3_raw[k,14], Bank3_raw[k,19]]
    Bank3_Sg_raw_5th[i,:]=[1080, 1155, 1790, 1775]*Bank3_raw[k,24]
end
# generate required data
Pd=transpose(hcat(0.95*Bank1_Pd_raw_5th/1000, 0.95*Bank2_Pd_raw_5th/1000,
    0.95*Bank3_Pd_raw_5th/1000,0.95*Bank1_Pd_raw_5th/1000, 0.95*Bank2_Pd_raw_5th/1000,
        0.95*Bank3_Pd_raw_5th/1000, 0.95*Bank1_Pd_raw_5th/1000, 0.95*Bank2_Pd_raw_5th/1000,
            0.95*Bank3_Pd_raw_5th/1000))
Sg_max=transpose(hcat(Bank1_Sg_raw_5th/1000, Bank2_Sg_raw_5th/1000,
    Bank3_Sg_raw_5th/1000, Bank1_Sg_raw_5th/1000, Bank2_Sg_raw_5th/1000,
        Bank3_Sg_raw_5th/1000, Bank1_Sg_raw_5th/1000, Bank2_Sg_raw_5th/1000,
            Bank3_Sg_raw_5th/1000))
###############################################################################
### Load Data from data directory
### Load Pd Data
# ### Load Price Data
Lambda_raw=CSV.read("data/Lambda_Feb1.csv");
Alpha_raw=CSV.read("data/Alpha_Feb1.csv");
Lambda=vec(Lambda_raw[:,1])
Alpha=vec(Alpha_raw[:,1])
Beta=15;

###############################################################################
### Performe the function
Pf=GML(Pd, Sg_max, Lambda, Beta, Alpha);

println("finished")
