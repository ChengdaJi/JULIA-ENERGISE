using JuMP
using Gurobi
using MAT
using LinearAlgebra
using CSV
using DataFrames

include("GML_struct.jl")
include("GML.jl")

###############################################################################
### Load Data
f_n=1;
e_n=288;

# demand_mat=matread("data/demand.mat");
# generation_mat=matread("data/generation.mat");
Price_mat=matread("data/price.mat");
# demand=demand_struct(demand_mat["Pd"][:,f_n:e_n],demand_mat["Qd"][:,f_n:e_n])
# Price=price_struct(Price_mat["alpha_10"][f_n:e_n],-15,
#     Price_mat["Lambda"][f_n:e_n]);
# Pg_trend=Pg_trend_struct(generation_mat["Pg1"][:,f_n:e_n],
#     generation_mat["Pg2"][:,f_n:e_n]);
# m=size(demand_mat["Pd"][:,f_n:e_n]);

# Pd_1=demand_mat["Pd"][:,f_n:e_n]/1000;
# println(size(Pd_1))
Pd_raw=CSV.read("Pd_Feb1.csv");
Pd=transpose(convert(Array, Pd_raw[:,:]))
# println(size(Pd))

# Qd=demand_mat["Qd"][:,f_n:e_n]/1000;
# Pg_trend_1=[generation_mat["Pg1"][:,f_n:e_n]; generation_mat["Pg2"][:,f_n:e_n]]
# println(size(Pg_trend_1))
Pg_raw=CSV.read("Pg_trend_Feb1.csv");
Pg_trend=transpose(convert(Array, Pg_raw[:,:]))
# Alpha=Price_mat["alpha_10"][f_n:e_n];
# Lambda=Price_mat["Lambda"][f_n:e_n];
Beta=15;
println(size(Pg_trend))

Lambda_raw=CSV.read("Lambda_Feb1.csv");
Alpha_raw=CSV.read("Alpha_Feb1.csv");
Lambda=vec(Lambda_raw[:,1])
Alpha=vec(Alpha_raw[:,1])
println(size(Alpha))
println(size(Lambda))
# println(length(Alpha))
# println(length(Lambda))
# println(dot(Alpha,Lambda))
#
GML(Pd, Pg_trend, Lambda, Beta, Alpha);
