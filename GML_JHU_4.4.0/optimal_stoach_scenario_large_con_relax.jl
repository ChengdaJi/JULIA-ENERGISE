function optimal_stoach_scenario_large_con_relax(current_time, obj, feedback, pd, pg,
    price, ancillary_type, base, network, T_emergency)
    println("===== GML - Optimization ")

    #################
    #  F,T or BN, T #
    #################
    ## The parameters

    T = obj.T; # number of time slots
    BN = obj.BN; # number of bus
    SN = obj.SN;
    Pg_min = obj.Pg_min; # minimum power generation

    Qf_max = obj.Qf_max; # maximum reactive power generation
    Qf_min = obj.Qf_min; # minimum eactive power generation

    B_max = obj.B_max; # maximum storage level
    B_min = obj.B_min; # minimum storage level

    R_max = obj.R_max; # the maximum discharge rate
    R_min = obj.R_min; # the minimum discharge rate

    delta_t = obj.delta_t; # time interval
    icdf = obj.icdf;

    # P_rsrv_max = obj.P_rsrv_max; # The maximum ancillary power

    P_rsrv_min = obj.P_rsrv_min; # The minimum ancillary power

    S=obj.S; # the apparents power

    V_min=obj.V_min;
    V_max=obj.V_max;

    beta=obj.beta; # price in cost function
    tau = obj.tau;

    r=obj.r; # the resistance
    x=obj.x; # the reactance
    k=obj.k; # time with ancellary

    Pd = pd.traj;

    B_feedback = feedback.B_feedback;
    P_rsrv_feedback = feedback.P_rsrv_feedback;

    bus = network.bus;
    branch = network.branch
    BrN=length(network.branch[:,1]);

    gen = reshape([1; 1; 51; 27; 0; 0; 1; 0.955; 10.67; 138; 1; 1.06; 0.96], 1,13);

    m = Model(with_optimizer(Mosek.Optimizer, QUIET=true, MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-3,
    MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-3))
    # define the real-time variables
    # Bus level
    @variable(m, Pg_rt[1:BN, 1])
    @variable(m, Qf_rt[1:BN, 1])
    @variable(m, B_rt[1:BN, 1])
    @variable(m, R_rt[1:BN, 1])
    @variable(m, P_bus_rt[1:BN, 1])
    @variable(m, Q_bus_rt[1:BN, 1])
    @variable(m, v_rt[1:BN, 1])

    # Branch level
    @variable(m, P_br_rt[1:BrN, 1])
    @variable(m, Q_br_rt[1:BrN, 1])
    @variable(m, l_br_rt[1:BrN, 1])
    # println(B_feedback)
    # println(" ---- Real Time Constraint Buildup ")
    for Bus=1:BN
        @constraint(m, 0<=Pg_rt[Bus,1]);
        @constraint(m, Pg_rt[Bus,1]<=
            positive_scalar(icdf*sqrt(pd.sigma[Bus,1]+pg.sigma[Bus,1])
            +pg.mu[Bus,1]));
        @constraint(m, R_min[Bus,1]<= R_rt[Bus,1]);
        @constraint(m, R_rt[Bus,1]<= R_max[Bus,1]);
        @constraint(m, B_rt[Bus,1]==B_feedback[Bus,1]);
        @constraint(m, Qf_min[Bus,1]<=Qf_rt[Bus,1]);
        @constraint(m, Qf_rt[Bus,1]<=Qf_max[Bus,1]);
        @constraint(m, V_min[Bus]^2<= v_rt[Bus,1]);
        @constraint(m, v_rt[Bus,1]<= V_max[Bus]^2);
        if iszero(findall(isagen->isagen==Bus, gen[:,1]))
            @constraint(m, P_bus_rt[Bus,1]==Pd[Bus,1]-Pg_rt[Bus,1]-R_rt[Bus,1]);
        else
            gen_branch_list=findall(isafbus->isafbus==Bus, branch[:,1]);
            @constraint(m,
            P_bus_rt[Bus,1]
            ==sum(P_br_rt[gen_branch_list[branch],1]
            for branch=1:length(gen_branch_list)));
        end

        if mod(Bus, 15) in [2,3,4]
            @constraint(m, [0.5*S, S, P_bus_rt[Bus,1],
                Q_bus_rt[Bus,1]] in RotatedSecondOrderCone());
        end
    end

    for Brh=1:BrN
        fbus = branch[Brh,1];
        tbus = branch[Brh,2];
        tbus_branch_list = findall(isafbus->isafbus==tbus, branch[:,1]);
        @constraint(m,P_br_rt[Brh,1]==
        r[Brh]*l_br_rt[Brh,1]+P_bus_rt[tbus, 1]
        +sum(P_br_rt[branch,1] for branch in tbus_branch_list));

        @constraint(m,Q_br_rt[Brh,1]==
        x[Brh]*l_br_rt[Brh,1]+Q_bus_rt[tbus, 1]
        +sum(Q_br_rt[branch,1] for branch in tbus_branch_list));

        @constraint(m,v_rt[fbus,1]==
        v_rt[tbus,1]+(r[Brh]^2+x[Brh]^2)*l_br_rt[Brh,1]
        -2*(r[Brh]*P_br_rt[Brh,1]+x[Brh]*Q_br_rt[Brh,1]));

        @constraint(m, [0.5*l_br_rt[Brh,1], v_rt[fbus,1],
        P_br_rt[Brh,1], Q_br_rt[Brh,1]] in RotatedSecondOrderCone());
    end

    @variable(m, Pg[1:SN, 1:BN, 1:T-1]);
    @variable(m, Qf[1:SN, 1:BN, 1:T-1]);
    @variable(m, B[1:SN, 1:BN, 1:T-1]);
    @variable(m, R[1:SN, 1:BN, 1:T-1]);
    @variable(m, P_bus[1:SN, 1:BN, 1:T-1])
    @variable(m, Q_bus[1:SN, 1:BN, 1:T-1])
    @variable(m, v[1:SN, 1:BN, 1:T-1])
    # Branch level
    @variable(m, P_br[1:SN,1:BrN, 1:T-1])
    @variable(m, Q_br[1:SN,1:BrN, 1:T-1])
    @variable(m, l_br[1:SN,1:BrN, 1:T-1])
    #
    for scenario = 1:SN
        for t=1:T-1
            for Bus=1:BN
                @constraint(m, Pg_min[Bus, t]<=Pg[scenario, Bus ,t]);
                @constraint(m, Pg[scenario, Bus,t]<=
                    positive_scalar(icdf*sqrt(pd.sigma[Bus,t+1]+pg.sigma[Bus,t+1])+pg.mu[Bus,t+1]));
                @constraint(m, B_min[Bus,t+1] <= B[scenario,Bus,t]);
                @constraint(m, B[scenario,Bus,t]<= B_max[Bus,t+1]);
                @constraint(m, R_min[Bus,t+1] <= R[scenario,Bus,t]);
                @constraint(m, R[scenario,Bus,t]<= R_max[Bus,t+1]);
                if t==1
                    @constraint(m, B[scenario,Bus,1] == B_rt[Bus,1]
                        -delta_t*(R_rt[Bus,1]*base.MVA))
                else
                    @constraint(m, B[scenario,Bus,t] ==
                        B[scenario,Bus,t-1] - R[scenario,Bus,t-1]*base.MVA*delta_t)
                end
                @constraint(m, Qf_min[Bus,t+1] <= Qf[scenario,Bus,t]);
                @constraint(m, Qf[scenario,Bus,t]<= Qf_max[Bus,t+1]);

                @constraint(m, V_min[Bus]^2<= v[scenario,Bus,t]);
                @constraint(m, v[scenario,Bus,t]<= V_max[Bus]^2);

                if iszero(findall(isagen->isagen==Bus, gen[:,1]))
                    @constraint(m, P_bus[scenario,Bus,t]==
                        Pd[Bus,t+1]-Pg[scenario,Bus,t]-R[scenario, Bus,t]);
                else
                    gen_branch_list=findall(isafbus->isafbus==Bus, branch[:,1]);
                    @constraint(m, P_bus[scenario,Bus,t]==
                    sum(P_br[scenario, gen_branch_list[branch],t]
                        for branch=1:length(gen_branch_list)));
                    # @constraint(m, P_bus[scenario,Bus,t]<=105/base.MVA)
                end
                if mod(Bus, 15) in [2, 3, 4]
                    @constraint(m, [0.5*S, S, P_bus[scenario,Bus,t],
                    Q_bus[scenario,Bus,t]] in RotatedSecondOrderCone());
                end
            end
            for Brh=1:BrN
                fbus = Int(branch[Brh,1]);
                tbus = Int(branch[Brh,2]);
                tbus_branch_list = findall(isafbus->isafbus==tbus, branch[:,1]);

                @constraint(m,P_br[scenario,Brh,t]==
                r[Brh]*l_br[scenario,Brh,t]+P_bus[scenario, tbus, t]+
                sum(P_br[scenario,branch,t] for branch in tbus_branch_list));

                @constraint(m,Q_br[scenario,Brh,1]==
                x[Brh]*l_br[scenario,Brh,t]+Q_bus[scenario,tbus, t]+
                sum(Q_br[scenario,branch,t]  for branch in tbus_branch_list));

                @constraint(m,v[scenario,tbus,t]==
                v[scenario,fbus,t]+(r[Brh]^2+x[Brh]^2)*l_br[scenario,Brh,t]-
                2*(r[Brh]*P_br[scenario,Brh,t]+x[Brh]*Q_br[scenario,Brh,t]));

                @constraint(m, [0.5*l_br[scenario,Brh,t], v[scenario,fbus,t],
                P_br[scenario,Brh,t], Q_br[scenario,Brh,t]] in RotatedSecondOrderCone());
            end
        end
    end

    ###########################################################################
    # variables and Constraints for reserve markets
    MVA=base.MVA
    if ancillary_type == "10min" || ancillary_type == "30min"
        Fake_Bus_list = [2]
        for bus = 3:BN
        # println(bus)
           if mod(bus, 15) in [2,3,4]
            # println("here")
               push!(Fake_Bus_list, bus)
           end
        end
        # println(Fake_Bus_list)
        Real_Bus_list = setdiff(2:BN, Fake_Bus_list)
        Non_affected_Bus_list = setdiff(Real_Bus_list, 5)
        # println(Real_Bus_list)
       # RT
        @variable(m, P_rsrv_rt)
        @variable(m, B_rsrv_rt)

        @constraint(m, P_rsrv_rt>=0)
        @constraint(m, B_rsrv_rt>=0)
       # @constraint(m, P_rsrv_rt==3)

        if current_time <= tau
            @constraint(m, B_rsrv_rt==0)

        elseif current_time - tau>=1
            ini_fb = max(current_time-tau-k+1,1);
            fni_fb = current_time-tau;
            length_fb = fni_fb - ini_fb +1;
            mult_fb = k-length_fb+1:k;
            # println(ini_fb)
            # println
            temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb]
            if length_fb >= 2
                for f_rsrv_fb_n=2:length_fb;
                    temp_f_rsrv_c_fb = (temp_f_rsrv_c_fb+
                        mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n]);
                end
            end
            @constraint(m, B_rsrv_rt==floor(delta_t*temp_f_rsrv_c_fb*1000)/1000)
            # println("nominal B_rsrv")
            # println(floor(delta_t*temp_f_rsrv_c_fb*1000)/1000)
            # println("sum B_feedback")
            # println(sum(B_feedback))
        end
        # for Bus = 2:BN
        #     list = setdiff(Real_Bus_list, Bus)
        #     # println(list)
        #     println(sum(B_feedback[bus, 1] for bus in list))
        #     @constraint(m, B_rsrv_rt <= sum(B_rt[bus, 1] for bus in list))
        # end

        # @constraint(m, B_rsrv_rt <= sum(B_rt))
        # println(sum(B_feedback[bus, 1] for bus in Non_affected_Bus_list))
        # println(B_feedback)
        # println(B_feedback)
        @constraint(m, B_rsrv_rt <=
            sum(B_rt[bus, 1] for bus in Non_affected_Bus_list))
        # @constraint(m, B_rsrv_rt <= 4)
        # Scenario
        @variable(m, P_rsrv[1:SN,1:T-1])
        @variable(m, B_rsrv[1:SN,1:T-1])
        # @constraint(m, P_rsrv==3)

        for scenario=1:SN
            for t_real=current_time+1:current_time+T-1
                @constraint(m, P_rsrv[scenario, t_real-current_time]>=0)
                @constraint(m, B_rsrv[scenario, t_real-current_time]>=0)
                ini = t_real-tau-k+1;
                fin = t_real-tau;
                if fin <=0
                    @constraint(m, B_rsrv[scenario, t_real-current_time]==0)
                elseif fin < current_time && fin>=1
                    ini_fb = max(1,ini);
                    fin_fb = fin;
                    length_fb = fin_fb-ini_fb+1;
                    mult_fb = k-length_fb+1:k;
                    temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb];
                    if length_fb >= 2
                        for f_rsrv_fb_n=2:length_fb;
                            temp_f_rsrv_c_fb = temp_f_rsrv_c_fb+
                                mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n];
                        end
                    end
                    @constraint(m, B_rsrv[scenario, t_real-current_time]==
                        floor(delta_t*temp_f_rsrv_c_fb*1000)/1000);
                elseif fin == current_time
                    if current_time == 1
                        @constraint(m, B_rsrv[scenario, t_real-current_time]==
                        delta_t*k*P_rsrv_rt);
                    else
                        ini_fb = max(1, ini);
                        fin_fb = fin-1;
                        length_fb = fin_fb-ini_fb+1;
                        mult_fb = k-length_fb:k-1;
                        temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb];
                        if length_fb >= 2
                            for f_rsrv_fb_n=2:length_fb;
                                temp_f_rsrv_c_fb = temp_f_rsrv_c_fb+
                                   mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n];
                             end
                        end
                        # @constraint(m, B_rsrv[scenario, t_real-current_time]==
                        #     floor(delta_t*k*P_rsrv_rt*1000)/1000);
                        @constraint(m, B_rsrv[scenario, t_real-current_time]==
                             delta_t*k*P_rsrv_rt
                             +floor(delta_t*temp_f_rsrv_c_fb*1000)/1000);
                    end
                 elseif ini < current_time && fin > current_time
                    if current_time == 1
                        ini_sc = current_time+1;
                        fin_sc = fin;
                        length_sc = fin_sc - ini_sc +1;
                        mult_sc = k-length_sc+1:k;
                        mult_rt = k-length_sc;
                         temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time]
                        if length_sc > 1
                            for f_rsrv_sc_n=2:length_sc
                                temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+
                                    mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
                            end
                        end
                        @constraint(m, B_rsrv[scenario, t_real-current_time] ==
                            delta_t*mult_rt*P_rsrv_rt
                            +delta_t*temp_f_rsrv_c_sc);
                    else
                        ini_sc = current_time+1;
                        fin_sc = fin;
                        length_sc = fin_sc - ini_sc +1;
                        mult_sc = k-length_sc+1:k;
                        temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time]
                        if length_sc > 1
                            for f_rsrv_sc_n=2:length_sc
                                temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+
                                    mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
                            end
                        end
                        mult_rt = k-length_sc;
                        ini_fb = max(1, ini);
                        fin_fb = current_time-1;
                        length_fb = fin_fb-ini_fb+1;
                        temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb];
                        if length_fb >= 2
                            for f_rsrv_fb_n=2:length_fb;
                                temp_f_rsrv_c_fb = temp_f_rsrv_c_fb+
                                    mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n];
                            end
                        end
                        @constraint(m, B_rsrv[scenario, t_real-current_time] ==
                        delta_t*mult_rt*P_rsrv_rt+
                        delta_t*temp_f_rsrv_c_sc+
                        floor(delta_t*temp_f_rsrv_c_fb*1000)/1000);
                    end
                elseif ini == current_time
                    ini_sc = current_time+1;
                    fin_sc = fin;
                    mult_sc = 2:k;
                    temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time];
                    for f_rsrv_sc_n=2:k-1
                        temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+
                            mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
                    end
                    @constraint(m, B_rsrv[scenario, t_real-current_time] ==
                        delta_t*P_rsrv_rt
                        +delta_t*temp_f_rsrv_c_sc);
                elseif ini > current_time
                    ini_sc = ini;
                    fin_sc = fin;
                    mult_sc = 1:k;
                    temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time]
                    for f_rsrv_sc_n=2:k
                        temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+
                            mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
                    end
                    @constraint(m, B_rsrv[scenario, t_real-current_time]
                        == delta_t*temp_f_rsrv_c_sc);

                    if  t_real-current_time >= T-tau-1
                        @constraint(m, P_rsrv[scenario, t_real-current_time]==0)
                    end
                end
                # for Bus = 2:BN
                #     list = setdiff(Real_Bus_list, Bus)
                #     @constraint(m, B_rsrv[scenario, t_real-current_time]
                #         <= sum(B[scenario, bus, t_real-current_time] for bus in list))
                # end
                if t_real in current_time+1:T_emergency+tau+k-1
                    @constraint(m, B_rsrv[scenario, t_real-current_time]
                        <= sum(B[scenario, bus, t_real-current_time] for
                        bus in Non_affected_Bus_list))
                else
                    for Bus in Non_affected_Bus_list
                        list = setdiff(Non_affected_Bus_list, Bus)
                        @constraint(m, B_rsrv[scenario, t_real-current_time]
                            <= sum(B[scenario, bus, t_real-current_time] for
                            bus in list))
                    end
                end
                # @constraint(m, B_rsrv[scenario, t_real-current_time]
                #      <= sum(B[scenario, :, t_real-current_time]))
                # @constraint(m, B_rsrv[scenario, t_real-current_time]
                #      <= 4)
           end
       end
    end
    if ancillary_type == "10min" || ancillary_type == "30min"
        @objective(m, Min,
            fn_cost_RHC_anc(delta_t,P_bus_rt,P_bus,Pg_rt,Pg,P_rsrv_rt,
            P_rsrv,price,pg,pd,beta,SN,obj,base,gen))
    # else
    #     @objective(m, Min,
    #         fn_cost_RHC_rt(delta_t,P_bus_rt,P_bus,Pg_rt,Pg,price,
    #             pg,pd,beta,SN,obj,base,gen))
    end

    status=optimize!(m);

    println(string("    ----", termination_status(m)))
    terminate_s = termination_status(m);
    # println(MOI.PrimalStatus())
    # println(MOI.DualStatus())
    cost_o = JuMP.objective_value(m);
    time_solve=MOI.get(m, MOI.SolveTime());
    println(string("    ----Solve Time: ", time_solve))
    println(string("    ----Optimal cost of whole horizon: ", cost_o))
    ## obtaining value

    ###############
    R_rt_o=JuMP.value.(R_rt)
    R_o=JuMP.value.(R)
    R_traj = hcat(R_rt_o, R_o[1,:,:])
    ################
    B_rt_o=JuMP.value.(B_rt)
    B_o=JuMP.value.(B)
    B_traj = hcat(B_rt_o, B_o[1,:,:])
    ###############
    Pg_rt_o=JuMP.value.(Pg_rt)
    Pg_o=JuMP.value.(Pg)
    Pg_traj = hcat(Pg_rt_o, Pg_o[1,:,:])

    ############
    P_bus_rt_o=JuMP.value.(P_bus_rt)
    P_bus_o=JuMP.value.(P_bus)
    P_bus_traj = hcat(P_bus_rt_o, P_bus_o[1,:,:])
    ############
    if ancillary_type == "10min" || ancillary_type == "30min"
        P_rsrv_rt_o=JuMP.value(P_rsrv_rt);
        P_rsrv_s=JuMP.value.(P_rsrv);
        P_rsrv_total = hcat(P_rsrv_rt_o[1,1], reshape(P_rsrv_s[1,:], 1, 287));
        B_rsrv_rt_o=JuMP.value(B_rsrv_rt);
        B_rsrv_s=JuMP.value.(B_rsrv);
        B_rsrv_total = hcat(B_rsrv_rt_o[1,1], reshape(B_rsrv_s[1,:], 1, 287));
    else
        P_rsrv_rt_o = 0;
        P_rsrv_total = zeros(1,T);
        B_rsrv_rt_o = 0;
    end

    ############
    Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
        +beta*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
        delta_t*P_rsrv_rt_o*price.alpha_ct
    # println(price.alpha_ct)

    alpha_1 = price.alpha_scenario[1,:]
    lambda_1 = price.lambda_scenario[1,:]

    println(string("    ----Optimzal cost at this instance: ", Cost_real))
    pg_upper = hcat(
    icdf*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])+pg.mu[:,:])

    return val_opt = (R=(R_rt_o), B=(B_rt_o), Pg=(Pg_rt_o),
        P_bus=(P_bus_rt_o), P_rsrv=(P_rsrv_rt_o), B_rsrv=(B_rsrv_rt_o),
        Cost_real=(Cost_real), time_solve=(time_solve),
        P_rsrv_total=(P_rsrv_total), B_rsrv_total=(B_rsrv_total), alpha_1=(alpha_1),
        R_traj = (R_traj), B_traj=(B_traj), Pg_traj=(Pg_traj), P_bus_traj=(P_bus_traj),
        P0_traj = (P_bus_traj[1,:]), pg_upper=(pg_upper), lambda_1=(lambda_1),
        pd=(pd), pg=(pg), terminate_s = (terminate_s))
end
