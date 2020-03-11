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

include("GML_large_emergency.jl")
include("traj_gen_br_mult_large.jl")
include("GML_struct.jl")
include("GML_RHC_Br_sto.jl")
include("GML_Networkbuildup.jl")
include("GML_RHC_Emergency.jl")
include("optimal_stoach_scenario_large_con_relax.jl")

# # of timeslots
T=288;
# # of secnarios
SN=6;


# ====================================================================
# M_list=[1 20];
# multiplier=1;
# ====================================================================
################################################################################

# Ancillar Markets Considered: "without", "10min", "30min"
ancillary_type = "10min"
# penetation levels
# p_rate = [0.25 0.5 0.75 1]; #[0.25, 0.5, 0.75, 1];
p_rate = 1;

# [90.95,99]
# icdf = [-1.2816, -1.6449, -2.3263]
icdf = -1.6449;
# prediction length [2hr, 1hr, .5hr]
pred_length = 24;

# max solar error
solar_error_max = 0.01;

# max solar error
B_cap = 3;

T_emergency = 3;
pcent_loss = 0.25;
################################################################################
multiplier = 1;
p_rate = 1;
# load price data
# for multiplier in M_list
	price_raw = read_price_data()
	delta_rt_raw=matread("../data/price_prediction.mat");

	# read demand data
	raw_data_mult = pd_pg_mult();

	# pd_raw = read_demand_data()
	pd_raw = read_NY_demand_data(raw_data_mult, multiplier)
	pd_noise = matread("../data/demand_noise.mat")["demand_noise"];


	# read pg noise
	pg_noise = matread("../data/solar_noise_0025.mat")["solar_noise"];

	# read pg_raw
	# pg_raw = read_solar_data()
	pg_raw=read_NY_solar_data(raw_data_mult, multiplier)


	baseMVA = 100;
	baseKV = 69;
	baseKA = baseMVA/baseKV;
	baseOhm = baseKV/baseKA;
	base=(MVA = (baseMVA), KV = (baseKV), KA = (baseKA), Ohm = (baseOhm));


	network = bbgm_generation(multiplier);
	BN = length(network.bus[:,1]);



val_opt = GML_large_emergency(ancillary_type, T, BN, SN,
	    p_rate, icdf, pred_length, solar_error_max, B_cap,
	    price_raw, delta_rt_raw, pd_raw, pd_noise, pg_noise, pg_raw,
	    base, multiplier, network, T_emergency, pcent_loss);
println("finish")
# pg=val_opt.pg
# pd=val_opt.pd

plot(1:288, reshape(val_opt.lambda_1,288,1), label="lambda1", linewidth=2)
# plot!(1:288, reshape(sum(pd.traj, dims=1)*100,288,1), label="Pd", linewidth=2)
# plot!(1:288, reshape(sum(pg.mu*100, dims=1),288,1), label="Pg aval", linewidth=2)
plot!(1:288, reshape(sum(val_opt.pg_upper, dims=1)*110,288,1), label="Pg upper", linewidth=2)
plot!(1:288, reshape(val_opt.P0_traj*100,288,1), label="P0", linewidth=2)
plot!(1:288, reshape(sum(val_opt.Pg_traj, dims=1)*100,288,1), label="Pg", linewidth=2)
# plot!(1:288, reshape(sum(val_opt.R_traj, dims=1)*100,288,1), label="R", linewidth=2)
# plot!(1:288, reshape(sum(val_opt.B_traj, dims=1),288,1), label="B", linewidth=2)
# #

# println(sum(val_opt.P_rsrv_total))


# plot(1:288, reshape(val_opt.alpha_1,288,1), label="alpha1", linewidth=2)
# plot!(1:288, reshape(val_opt.B_rsrv_total,288,1), label="B_rsrv", linewidth=2)
# plot!(1:288, reshape(val_opt.P_rsrv_total,288,1), label="P_rsrv", linewidth=2)
# plot!(1:288, reshape(val_opt.lambda_1,288,1), label="lambda1", linewidth=2)
