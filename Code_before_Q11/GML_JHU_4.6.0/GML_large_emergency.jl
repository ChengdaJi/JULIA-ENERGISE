function GML_large_emergency(ancillary_type, T, BN, SN,
    p_rate, icdf, pred_length, solar_error_max, B_cap,
    price_raw, delta_rt_raw, pd_raw, pd_noise, pg_noise, pg_raw,
    base, multiplier, network, T_emergency, pcent_increse)
    k=12;
    for current_time=T_emergency:T
    # current_time=1
        ct_printout = string("===== GML - At Time ", current_time);
        println("=================================================")
        println(ct_printout)
        if current_time == T_emergency
            ################################
            # old_feedback = read_old_data(T_emergency,multiplier)
            # B_feedback = old_feedback.B_feedback
            # P_rsrv_feedback = old_feedback.P_rsrv_feedback
            ################################
            # P_rsrv_feedback = [];
            P_rsrv_feedback = hcat(zeros(1, T_emergency-3), 0.1*ones(1,2));
            B_feedback_one = vcat(zeros(3,1), B_cap/12*ones(12,1));
            B_feedback = vcat(0, repeat(B_feedback_one, multiplier));

        else
            B_feedback = read_B_out()
            P_rsrv_feedback = read_RSRV_out()
        end

        price = price_traj(current_time, ancillary_type,
            price_raw, delta_rt_raw, T, pred_length);

        pd = pd_traj_large(current_time, pd_raw, pd_noise, T, pred_length, base,multiplier);

        feedback = (B_feedback=(B_feedback), P_rsrv_feedback=(P_rsrv_feedback));

        pg = pg_traj_large(current_time, pg_raw, pg_noise,
            solar_error_max, p_rate, T, pred_length,  base, multiplier);
        println(sum(pg.mu))

        obj = GML_Sys_Ava_large(T, BN, SN, pd, ancillary_type, icdf,
            B_cap, base, network)

        if current_time in T_emergency:T_emergency+k-1

            # val_opt = optimal_stoach_scenario_large_con(current_time, obj, feedback,
            #     pd, pg, price, ancillary_type, base, network);
            val_opt = optimal_stoach_scenario_large_con_solar(current_time, obj, feedback,
                pd, pg, price, ancillary_type, base, network, pcent_increse);
        else
            val_opt = optimal_stoach_scenario_large_con(current_time, obj, feedback,
                pd, pg, price, ancillary_type, base, network);
        end

        IF_OP=string(val_opt.terminate_s)
        ############################
        # if IF_OP == "INFEASIBLE"
        #     println(">>>>>>>>>>>>>>>>>>>")
        #     break
        # end
        # ###########################
        B_feedback_out=zeros(BN,1)
        for bus = 1:BN
            B_feedback_out[bus, 1] = val_opt.B[bus,1] -
                floor(val_opt.R[bus,1]/12*base.MVA*1000)/1000;
        end
        # println(val_opt.P_rsrv)
        if val_opt.P_rsrv <= 0.001
            P_rsrv_feed = 0.0;
        else
            P_rsrv_feed = floor(val_opt.P_rsrv*1000)/1000
        end

        if current_time == 1
            P_rsrv_feedback_temp = [P_rsrv_feed]
        else
            # println(size(P_rsrv_feedback))
            P_rsrv_feedback_temp = hcat(P_rsrv_feedback, P_rsrv_feed)
            # push!(P_rsrv_feedback,val_opt.P_rsrv)
        end
        # println(B_feedback_out)
        # println(P_rsrv_feedback_temp)
        write_B_out(B_feedback_out)
        write_RSRV_out(P_rsrv_feedback_temp)


        mkpath("./result")
        write_output_out(val_opt,P_rsrv_feed,
            string("./result/M", multiplier, "P_rate", p_rate, "Time", current_time, ".csv"))
    end
    # return val_opt
end

# plot(1:288, reshape(val_opt.lambda1,288,1), label="lambda1", linewidth=2)
# plot!(1:288, reshape(sum(pd.traj, dims=1)*100,288,1), label="Pd", linewidth=2)
# # plot!(1:288, reshape(sum(pg.mu, dims=1)*100,288,1), label="Pg aval", linewidth=2)
# # plot!(1:288, reshape(val_opt.pg_upper*110,288,1), label="Pg upper", linewidth=2)
# plot!(1:288, reshape(val_opt.P0[1,:]*100,288,1), label="P0", linewidth=2)
# # plot!(1:288, reshape(sum(val_opt.Pg, dims=1)*100,288,1), label="Pg", linewidth=2)
# # plot!(1:288, reshape(sum(val_opt.R, dims=1)*100,288,1), label="R", linewidth=2)
# plot!(1:288, reshape(sum(val_opt.B, dims=1),288,1), label="B", linewidth=2)
