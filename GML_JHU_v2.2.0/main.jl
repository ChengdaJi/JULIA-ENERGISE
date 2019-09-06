# this version is for trade-off analysis
#

using MAT
using LinearAlgebra
using JuMP
using Mosek
using MosekTools
using CSV
using DataFrames
using Plots

include("traj_gen.jl")
include("GML_struct.jl")
include("GML_RHC.jl")
## sys para
ancillary_type="10min"
# ancillary_type="30min"
# ancillary_type="Without"
# # of timeslots
T=288;
# # of banks
# BN=252;
BN=3;
# # of feeders
F=BN*4;
# # of secnarios
SN=6;

################################################################################
# penetation level
p_rate=0.5; # 0.25 0.5 0.75 1

icdf = -1.6449; #0.95
# icdf = -1.2816 #0.9
# icdf = -1.0364 #0.85

# # prediction Pred_length
Pred_length=24; # 5 min per slots
## read solar
solar_error_max = 0.025;
################################################################################

# # initial SOC
# load price data
price_raw = read_price_data()
delta_rt_raw=matread("../data/price_prediction.mat");

## read read_demand_data
pd_raw = read_demand_data()
pd_noise = matread("../data/demand_noise.mat")["demand_noise"];



if solar_error_max == 0.025
    pg_noise = matread("../data/solar_noise_0025.mat")["solar_noise"];
elseif solar_error_max == 0.05
    pg_noise = matread("../data/solar_noise_005.mat")["solar_noise"];
elseif solar_error_max == 0.1
    pg_noise = matread("../data/solar_noise_01.mat")["solar_noise"];
end
pg_raw = read_solar_data()
# #
for current_time=1:T
    ct_printout = string("===== GML - At Time ", current_time);
    println("=================================================")
    println(ct_printout)
    if current_time ==1
        global P_rsrv_feedback = [];
        global B_feedback=reshape([0.0910,0.0942,0.1434,0.1402,0.1382,0.2144,0.1192,0.1813,
            0.1355,0.0664,0.1182,0.0580],12,1);
    end
    global feedback = (B_feedback=(B_feedback), P_rsrv_feedback=(P_rsrv_feedback));
    price = price_traj(current_time, ancillary_type, price_raw, delta_rt_raw, T, Pred_length);
    # plot(1:287, price.alpha_scenario[1,:])
    pd = pd_traj(current_time, pd_raw, pd_noise, BN, T, Pred_length);
    pg = pg_traj(current_time, pg_raw, pg_noise, solar_error_max, p_rate, T, Pred_length);
    obj = GML_Sys_Ava(T, F, BN, SN, pd, ancillary_type, icdf);
    val_opt = optimal_stoach_scenario(current_time, obj, feedback, pd, pg,
    price, ancillary_type);
    for feeder = 1:12
        B_feedback[feeder, 1] = val_opt.B[feeder,1] - val_opt.R[feeder,1]/12;
    end
    if current_time == 1
        P_rsrv_feedback = [val_opt.P_rsrv]
    else
        push!(P_rsrv_feedback,val_opt.P_rsrv)
    end
    write_output_out(val_opt, current_time)
    println("=================================================")
end



# Pd_rt=pd_raw_one["Pd"][:,1:288]/1000 + error_fixed;
# println("GML - optimization finished")
