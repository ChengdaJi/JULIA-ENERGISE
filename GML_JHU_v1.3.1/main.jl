# this version is for FOL-GML integration

using MAT
using LinearAlgebra
using Statistics
using JuMP
using Mosek
using MosekTools
using CSV
using DataFrames
using Plots

include("traj_gen_det.jl")
include("GML_struct.jl")
include("GML_RHC.jl")
## sys para
# current time slot 1~288
current_time = 1
# Battery SOC
B_input = 0.1;
# run the GML_RHC_DEG function
GML_setpoint = GML_RHC_DET(current_time, B_input)
# Optimal battery discharge
GML_setpoint.R
println(string("Optimal battery discharge value ", GML_setpoint.R, " MW"))
# Optimal Solar Generation
GML_setpoint.Pg
println(string("Optimal solar generation ", GML_setpoint.Pg, " MW"))
# Optimal Agg Reactive Power
GML_setpoint.Qf
println(string("Optimal aggergated reactive power ", GML_setpoint.Qf, " MVA"))
# optimal headnode voltage
GML_setpoint.v
println(string("Optimal headnode voltage", GML_setpoint.v, " KV"))
