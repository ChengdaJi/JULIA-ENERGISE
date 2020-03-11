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

function GML(ancillary_type, T, BN, F, SN,
    p_rate, icdf, pred_length, solar_error_max, B_cap,
    price_raw, delta_rt_raw, pd_raw, pd_noise, pg_noise, pg_raw,
    folder, filename)

    for current_time=1:1
        ct_printout = string("===== GML - At Time ", current_time);
        println("=================================================")
        println(ct_printout)
        if current_time == 1
            global P_rsrv_feedback = [];
            global B_feedback=reshape([0.0910,0.0942,0.1434,0.1402,0.1382,0.2144,0.1192,0.1813,
                0.1355,0.0664,0.1182,0.0580],12,1);
        end

        feedback = (B_feedback=(B_feedback), P_rsrv_feedback=(P_rsrv_feedback));
        price = price_traj(current_time, ancillary_type, price_raw, delta_rt_raw, T, pred_length);
        # println(size(price.alpha_scenario))
        # println(size(price.lambda_scenario))
        pd = pd_traj(current_time, pd_raw, pd_noise, BN, T, pred_length);
        # println(size(pd.traj))
        # println(size(pd.sigma))
        pg = pg_traj(current_time, pg_raw, pg_noise, solar_error_max, p_rate, T, pred_length);
        # println(size(pg.mu))
        # println(size(pg.sigma))
        obj = GML_Sys_Ava(T, F, BN, SN, pd, ancillary_type, icdf, B_cap);
        val_opt = optimal_stoach_scenario(current_time, obj, feedback, pd, pg,
            price, ancillary_type);
        # for feeder = 1:12
        #     B_feedback[feeder, 1] = val_opt.B[feeder,1] - val_opt.R[feeder,1]/12;
        # end
        # if current_time == 1
        #     P_rsrv_feedback = [val_opt.P_rsrv]
        # else
        #     push!(P_rsrv_feedback,val_opt.P_rsrv)
        # end
        # mkpath(folder)
        # write_output_out(val_opt,
        #     string(folder, filename, "_time", current_time, ".csv"))
        # println("=================================================")

    end
end
