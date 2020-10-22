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
date = "aug01"

#############################################################################
raw_data = pre_process(date)

##############################################################################
initial_time = 12*12+1;
end_time = 13*12;
for current_time = initial_time:end_time

	##################
	if current_time == initial_time
		B_inp = [0.2935, 0.4566, 0.80902]/2;
	else
		B_inp = [0; 0; 0];
	end
	# 3*1 matrix, where B_inp[1] represent the B soc at feeder 1.
	# Unit MWhr

	##################

	solar_percent = 1;
	battery_capacity_percent = 1;

	val_opt = GML_solver(date, current_time, B_inp, raw_data, solar_percent, battery_capacity_percent)
	write_R_FOL_out(current_time, val_opt)
	write_Pg_FOL_out(current_time, val_opt)
	write_Q_FOL_out(current_time, val_opt)
	write_P_FOL_out(current_time, val_opt)
	write_v_FOL_out(current_time, val_opt)
#

end