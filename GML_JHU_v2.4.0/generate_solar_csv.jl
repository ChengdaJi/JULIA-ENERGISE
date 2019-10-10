using CSV
using DataFrames
using MAT
using Plots

include("traj_gen_12_feeders.jl")
include("GML_struct.jl")
include("GML_RHC.jl")
include("ancillary_functions.jl")

p_rate=1;
BN=6;
T=288;
Pred_length=12;
pd_raw = read_demand_data()
pg_raw_one=matread("../data/solar.mat");
pg = pg_traj_raw(1, pg_raw_one, pd_raw.pd_rt[:,1:288], p_rate, BN,T, Pred_length);
solar_output(pg.mu)
pg_raw = read_solar_data()
# plot(1:288, reshape())
