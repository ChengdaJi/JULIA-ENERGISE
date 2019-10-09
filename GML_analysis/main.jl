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
using Statistics

include("traj_gen_det.jl")
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
SN=1;

################################################################################
# penetation level
Solar_pene=[0.25 0.5 0.75 1]; # 0.25 0.5 0.75 1

# deterministic
icdf = 0;
Battery_sizes = [3 15 30];
#
################################################################################

# # initial SOC
# load price data
price_raw = read_price_data()
delta_rt_raw=matread("../data/price_prediction.mat");

## read read_demand_data
pd_raw = read_demand_data()
pd_noise = matread("../data/demand_noise.mat")["demand_noise"];



# if solar_error_max == 0.025
#     pg_noise = matread("../data/solar_noise_0025.mat")["solar_noise"];
# elseif solar_error_max == 0.05
#     pg_noise = matread("../data/solar_noise_005.mat")["solar_noise"];
# elseif solar_error_max == 0.1
#     pg_noise = matread("../data/solar_noise_01.mat")["solar_noise"];
# end
pg_raw = read_solar_data()
# #
current_time=1;
for p_rate in Solar_pene
    for B_input in Battery_sizes
        ct_printout = string("===== Solar ", p_rate, " Battery ", B_input);
        println("=================================================")
        println(ct_printout)
        P_rsrv_feedback = [];
        B_feedback=reshape([0.0910,0.0942,0.1434,0.1402,0.1382,0.2144,0.1192,0.1813,
                0.1355,0.0664,0.1182,0.0580],12,1);
        feedback = (B_feedback=(B_feedback), P_rsrv_feedback=(P_rsrv_feedback));
        price = price_traj_det(current_time, ancillary_type, price_raw, T);
        pd = pd_traj_det(current_time, pd_raw, T)
        pg = pg_traj_det(current_time, pg_raw, p_rate, T);
        # println(var((pd.da-pd.traj)))
        # p
        # temp=((pg.da-pg.mu)./pg.mu)
        # println(var(temp))

        obj = GML_Sys_Ava(T, F, BN, SN, pd, ancillary_type, icdf, B_input);
        val_opt = optimal_stoach_scenario(current_time, obj, feedback, pd, pg,
        price, ancillary_type);
        pg_aval = sum(pg.mu, dims=1)
        println(string("Solar Utlizing rate ", sum(val_opt.Pg)/sum(pg_aval)));
        println(string("Total avaliability is ", sum(pg_aval)));
        println(string("Total demand is ", sum(pd.traj)));
        println(string("Total solar is ", sum(val_opt.Pg)));
        # write_output_out(val_opt, current_time)
        println("=================================================")
    end
end



# Pd_rt=pd_raw_one["Pd"][:,1:288]/1000 + error_fixed;
# println("GML - optimization finished")
