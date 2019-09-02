using MAT
using LinearAlgebra
using JuMP
using Mosek
using MosekTools

include("traj_gen_12_feeders.jl")
include("GML_struct.jl")
include("GML_RHC.jl")
## sys para
ancillary_type="without"
# # of timeslots
T=288;
# # of banks
# BN=252;
BN=3;
# # of feeders
F=BN*4;
# # of secnarios
SN=6;
# # initial SOC
B_feedback_1=reshape([0.0910,0.0942,0.1434,0.1402,0.1382,0.2144,0.1192,0.1813,
    0.1355,0.0664,0.1182,0.0580],12,1);
# B_feedback_c=Array(Array{Float64,2}, BN/3);
B_feedback_c=[];
for i=1:BN/3
    push!(B_feedback_c, B_feedback_1)
end
B_feedback=vcat(B_feedback_c...)
P_rsrv_feedback = [];

pd_raw_one = matread("../data/demand.mat")
pd_MAPE = matread("../data/ForecastError.mat")
pg_raw_one=matread("../data/solar.mat");
pg_RMSE = matread("../data/ForecastError.mat");
price_raw=matread("../data/price_data.mat");
delta_rt_raw=matread("../data/price_prediction.mat");
error_fixed=matread("../data/error.mat")["e"];
p_rate=0.5;

for current_time=1:1
    price = price_traj(current_time, ancillary_type, price_raw, delta_rt_raw);
    println(size(price.probability))
    pd=pd_traj(current_time, pd_raw_one, pd_MAPE, error_fixed[:,current_time],BN);
    println(size(pd.traj))
    println(size(pd.sigma))
    pg=pg_traj(current_time, pg_raw_one, pg_RMSE, pd.traj,p_rate,BN);
    println(size(pg.mu))
    println(size(pg.mu_scenario))
    println(size(pg.mu_rt))
    println(size(pd.sigma))
    B_rate=3;
    R_rate=1/3;
    obj = GML_Sys_Ava(T, F, BN, pd, B_rate, R_rate, ancillary_type);
    optimal_stoach_scenario_offline(current_time, obj, pd, pg,
    price, SN, B_feedback ,P_rsrv_feedback, ancillary_type);
end



Pd_rt=pd_raw_one["Pd"][:,1:288]/1000 + error_fixed;
println("GML - optimization finished")
