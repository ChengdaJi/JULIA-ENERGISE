# this version is for trade-off analysis
#

using MAT
using LinearAlgebra
using JuMP
using Mosek
using MosekTools
using CSV
using DataFrames
# using Plots
using Statistics
using Ipopt
using PowerModels

# only for ploting purpose, no need to use
# using PyCall
# pygui(true)
# using PyPlot

include("GML_NYISO.jl")
include("demand_solar_price.jl")
include("NYISO_read.jl")

# data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
# data_trace_shunt=CSV.File("../data/NYISO-data/shunt_ORU.csv") |> DataFrame
# # data_trace_branch=CSV.File("../data/NYISO-data/branch.csv") |> DataFrame
# data_trace_branch=CSV.File("../data/NYISO-data/branch_edit.csv") |> DataFrame
# data_trace_gen=CSV.File("../data/NYISO-data/gen.csv") |> DataFrame


###################
# specify the date of interest aug01-aug05;
date = "aug02"

#############################################################################
raw_data = pre_process(date)

##############################################################################
current_time = 288;

##################
B_inp = [0; 0; 0];
# 3*1 matrix, where B_inp[1] represent the B soc at feeder 1.
# Unit MWhr

val_opt = GML_solver(date, current_time, B_inp, raw_data)
#
println("battery discharge [p.u.] size 3 feeder * 12 time slots (1 hr)")
println(val_opt.R_FOL)
println("Solar usage [p.u.] size 3 feeder * 12 time slots (1 hr)")
println(val_opt.Pg_FOL)
println("Net real power (withdraw) [p.u.] size 3 feeder * 12 time slots (1 hr)")
println(val_opt.P_FOL)
println("Net reactive power (withdraw) [p.u.] size 3 feeder * 12 time slots (1 hr)")
println(val_opt.Q_FOL)
println("voltage magnitude [p.u.] size 3 feeder * 12 time slots (1 hr)")
println(val_opt.v_FOL)
#
# # write_P_gen_out(val_opt.P_gen_traj)
# # write_P_bus_out(val_opt.P_bus_traj)
#
# println("finish")


# folder = string("t=", current_time)
# mkdir(folder)
# write_branch_real_output(current_time,val_opt)
# write_branch_reactive_output(current_time,val_opt)
# write_bus_real_output(current_time,val_opt)
# write_bus_reactive_output(current_time,val_opt)
# write_bus_voltage_output(current_time,val_opt)
# write_generator_real_output(current_time,val_opt)
# write_generator_reactive_output(current_time,val_opt)
# # # #
# # println("=================================================")
# # # # plot(1:288, reshape(val_opt.lambda_1,288,1), label="lambda1", linewidth=2)
# # # # plot!(1:288, reshape(sum(pd.traj, dims=1)*100,288,1), label="Pd", linewidth=2)
# # # # plot!(1:288, reshape(sum(pg.mu, dims=1)*100,288,1), label="Pg aval", linewidth=2)
# # # # plot!(1:288, reshape(sum(val_opt.pg_upper, dims=1)*110,288,1), label="Pg upper", linewidth=2)
# # # # plot!(1:288, reshape(sum(val_opt.P0_traj, dims=1)*100,288,1), label="P0", linewidth=2)
# # Pg=reshape(sum(val_opt.Pg_traj, dims=1)*100,288,1);
# # Qg=reshape(sum(val_opt.Qg_traj, dims=1)*100,288,1);
# # plot(1:288, sqrt.(Pg.^2+Qg.^2))
# # plot!(1:288, reshape(sum(val_opt.Pg_traj, dims=1)*100,288,1), label="Pg", linewidth=2)
# # plot!(1:288, reshape(sum(val_opt.Qg_traj, dims=1)*100,288,1), label="Qg", linewidth=2)
# # # # plot!(1:288, reshape(sum(val_opt.R_traj, dims=1)*100,288,1), label="R", linewidth=2)
# # # # plot!(1:288, reshape(sum(val_opt.B_traj, dims=1),288,1), label="B", linewidth=2)
# # #
# # # # plot(1:576, reshape(sum(pd_raw.pd_rt, dims=1),576,1), label="rt", linewidth=2)
# # # # plot(1:576, reshape(sum(pd_raw.pd_da, dims=1),576,1), label="da", linewidth=2)
# # # # plot(1:288, reshape(val_opt.P_rsrv_total,288,1), label="P_rsrv", linewidth=2)
# # # # plot!(1:288, reshape(val_opt.alpha_1,288,1), label="alpha", linewidth=2)
