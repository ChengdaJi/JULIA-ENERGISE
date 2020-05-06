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

include("GML_NYISO.jl")
include("traj_gen_NYISO.jl")
include("GML_struct.jl")
include("NYISO_read.jl")
# data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
# data_trace_shunt=CSV.File("../data/NYISO-data/shunt.csv") |> DataFrame
# data_trace_branch=CSV.File("../data/NYISO-data/branch.csv") |> DataFrame
# data_trace_gen=CSV.File("../data/NYISO-data/gen.csv") |> DataFrame

# sys para
# ancillary_type="10min"
# ancillary_type="30min"
ancillary_type="10min";
# # of timeslots
T=288;
# # of bus
# NoShunt=209;
# # of secnarios
SN=1;

baseMVA = 1;
################################################################################
# penetation level
# Solar_pene=[0]; # 0.25 0.5 0.75 1
# Solar_pene=[0.25]
# deterministic
icdf = 0;
# Battery_sizes = [3 15 30];
# Battery_sizes = [3];
#
solar_error_max = 0;
################################################################################
# # shunt_struct = shunt_data(data_trace_bus,data_trace_shunt);

bus_struct = bus_data();
branch_struct =branch_data();
gen_struct = gen_data();
NoBus = length(bus_struct.baseKV)
# # initial SOC
# load price data
price_raw = read_price_data()
delta_rt_raw=matread("../data/price_prediction.mat");


## read read_demand_data
pd_raw = read_NY_demand_data(bus_struct)
# println(sum(pd_raw.pd_rt, dims=1))
# println(sum(pd_raw.pd_rt, dims=2))
# println(sum(pd_raw.pd_rt))


pd_noise = matread("../data/demand_noise.mat")["demand_noise"];
#
if solar_error_max == 0.025
    pg_noise = matread("../data/solar_noise_0025.mat")["solar_noise"];
elseif solar_error_max == 0.05
    pg_noise = matread("../data/solar_noise_005.mat")["solar_noise"];
elseif solar_error_max == 0.1
    pg_noise = matread("../data/solar_noise_01.mat")["solar_noise"];
elseif solar_error_max == 0
    pg_noise = zeros(1,12);
end
# println("here")
pg_raw = read_NY_solar_data(bus_struct)

# println(size(pg_raw.pg_rt))
# println(sum(pg_raw.pg_rt))
# println(size(pg_raw.pg_da))
# println(sum(pg_raw.pg_da))
# plot(1:576, reshape(sum(pd_raw.pd_rt, dims=1),576,1), label="pd rt", linewidth=2)
# plot!(1:576, reshape(sum(pd_raw.pd_da, dims=1),576,1), label="pd da", linewidth=2)
# plot!(1:576, reshape(sum(pg_raw.pg_rt, dims=1),576,1), label="pg rt", linewidth=2)
# plot!(1:576, reshape(sum(pg_raw.pg_da, dims=1),576,1), label="pg da", linewidth=2)

multiplier = demand_multiplier(bus_struct, pd_raw)

###########################################################################
# #
current_time=10;
p_rate = 0;
# B_cap = 150;
B_cap = 0;
ct_printout = string("===== Solar ", p_rate, " Battery ", B_cap);
println("=================================================")
println(ct_printout)
if current_time == 1
    P_rsrv_feedback = []
else
    P_rsrv_feedback = zeros(current_time-1,1);
end
B_feedback=zeros(NoBus, 1);
feedback = (B_feedback=(B_feedback), P_rsrv_feedback=(P_rsrv_feedback));
price = price_traj_det(current_time, ancillary_type, price_raw, T);
################################################################################
# plotting prices
################################################################################
# # plot(1:288, reshape(price.lambda_scenario[1,:], 288,1), label="lambda")
# # plot!(1:288, reshape(price.alpha_scenario[1,:], 288,1), label="alpha")
# # lambda = reshape(price.lambda_scenario[1,:], 288,1)
# # alpha = reshape(price.alpha_scenario[1,:], 288,1)
# # price = hcat(lambda, alpha)
# # CSV.write("price.csv", DataFrame(price, [:lambda, :alpha]))
# # println(sum(price.alpha_scenario[1,:]))
################################################################################
pd = pd_traj_pu_det(current_time, pd_raw, T,  NoBus, baseMVA, multiplier);
################################################################################
# for ploting demands
################################################################################
# println(sum(pd.traj, dims=2))
# println(sum(pd.traj))
# println(size(pd.traj))
# println(size(pd.sigma))
# println(size(pd.ct))
# plot(1:288, reshape(sum(pd.traj, dims=1), 288,1), label="traj")
# plot!(1:288, reshape(sum(pd.sigma, dims=1), 288,1), label="sigma")
# println(maximum(pd.traj, dims=2))
################################################################################
pg = pg_traj_pu_det(current_time, pg_raw, p_rate, T, NoBus, baseMVA, multiplier);
################################################################################
# for ploting solars
################################################################################
# println(size(pg.mu))
# println(size(pg.sigma))
# println(size(pg.mu_ct))
# # println(sum(pg.mu)/sum(pd.traj))
# # plot(1:288, reshape(sum(pg.mu, dims=1), 288,1), label="traj")
################################################################################
obj = GML_Sys_Ava_NYISO(T, pd, ancillary_type, B_cap, icdf, bus_struct,
	branch_struct, gen_struct, baseMVA);

# val_opt = optimal_NYISO_SOCP(SN, current_time, obj, ancillary_type, baseMVA,
#     feedback, pd, pg, price, bus_struct, branch_struct, gen_struct);

val_opt = optimal_NYISO(SN, current_time, obj, ancillary_type, baseMVA,
    feedback, pd, pg, price, bus_struct, branch_struct, gen_struct);

write_branch_real_output(val_opt)
write_branch_reactive_output(val_opt)
write_bus_real_output(val_opt)
write_bus_reactive_output(val_opt)
write_bus_voltage_output(val_opt)
write_generator_real_output(val_opt)
write_generator_reactive_output(val_opt)

println("=================================================")
# plot(1:288, reshape(val_opt.lambda_1,288,1), label="lambda1", linewidth=2)
# plot!(1:288, reshape(sum(pd.traj, dims=1)*100,288,1), label="Pd", linewidth=2)
# plot!(1:288, reshape(sum(pg.mu, dims=1)*100,288,1), label="Pg aval", linewidth=2)
# plot!(1:288, reshape(sum(val_opt.pg_upper, dims=1)*110,288,1), label="Pg upper", linewidth=2)
# plot!(1:288, reshape(sum(val_opt.P0_traj, dims=1)*100,288,1), label="P0", linewidth=2)
# plot!(1:288, reshape(sum(val_opt.Pg_traj, dims=1)*100,288,1), label="Pg", linewidth=2)
# plot!(1:288, reshape(sum(val_opt.R_traj, dims=1)*100,288,1), label="R", linewidth=2)
# plot!(1:288, reshape(sum(val_opt.B_traj, dims=1),288,1), label="B", linewidth=2)

# plot(1:576, reshape(sum(pd_raw.pd_rt, dims=1),576,1), label="rt", linewidth=2)
# plot(1:576, reshape(sum(pd_raw.pd_da, dims=1),576,1), label="da", linewidth=2)
# plot(1:288, reshape(val_opt.P_rsrv_total,288,1), label="P_rsrv", linewidth=2)
# plot!(1:288, reshape(val_opt.alpha_1,288,1), label="alpha", linewidth=2)
