function GML_Sys_Ava_large(T, BN, SN, pd, ancillary_type, icdf, B_cap, base, network)


    mult = Int((BN-1)/15);
    println("===== GML - Boundaries Buildup");
    ###############################################################################
    # feeder level
    Qf_max=0.05*positive_array(pd.traj);
    Qf_min = -Qf_max;
    # minimum solar
    Pg_min = zeros(BN, T);
    # battery
    B_rate=B_cap;
    R_rate=1/3;
    B_max = vcat(0,
        repeat(vcat(zeros(3,T), ones(12,1)*B_rate*ones(1,T)/12), mult));

    B_min = zeros(BN,T);
    R_max = R_rate*B_max/base.MVA;
    R_min = -R_max;

    W=zeros(BN,T);
    # ancillary
    if ancillary_type == "10min"
        tau=2;
    elseif ancillary_type == "30min"
        tau=6;
    else
        tau=2;
    end
    P_rsrv_min=zeros(1,T);
    k=12;

    delta_t = 1/12;
    S=35/base.MVA;
    # V
    # Base_V=69;
    # V_min = 0.96;
    # V_max = 1.06;
    bus = network.bus
    branch = network.branch
    V_max = bus[:,12];
    # println(V_max)
    V_min = bus[:,13];
    # impedance
    r=branch[:,3];
    x=branch[:,4];
    # println(r)

    # penalty
    beta=-15;
    println("===== GML - finish building boundaries =====")
    obj=obj_struct(T,0,BN,SN, Pg_min,Qf_max,Qf_min,B_max,B_min,R_max,R_min,W,
    delta_t,P_rsrv_min,tau,k,V_max,V_min,r,x,beta,S,icdf);
    return obj
end


function optimal_stoach_scenario_large(current_time, obj, feedback, pd, pg, price, ancillary_type, base, network)
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

        @constraint(m, B_rsrv_rt <= sum(B_rt))
        # println(sum(B_feedback[bus, 1] for bus in Non_affected_Bus_list))
        # @constraint(m, B_rsrv_rt <=
        #     sum(B_rt[bus, 1] for bus in Non_affected_Bus_list))
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


                @constraint(m, B_rsrv[scenario, t_real-current_time]
                     <= sum(B[scenario, :, t_real-current_time]))
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
    # println
    if sum(Pg_rt_o)<=sum(pg.mu_ct)
    ############
        Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
            +beta*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    else
        Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
            +price.lambda_ct*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    end
    # println(price.alpha_ct)
    P_cul = sum(pg.mu_ct)-sum(Pg_rt_o)
    println("curtailment solar")
    println(P_cul)
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
        pd=(pd), pg=(pg), terminate_s = (terminate_s), P_cul = (P_cul))
end

function optimal_stoach_scenario_large_con(current_time, obj, feedback, pd, pg, price, ancillary_type, base, network)
    println("===== GML - Optimization with N-1")

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
        for Bus in Real_Bus_list
            list = setdiff(Real_Bus_list, Bus)
            # println(list)
            # println(sum(B_feedback[bus, 1] for bus in list))
            @constraint(m, B_rsrv_rt <= sum(B_rt[bus, 1] for bus in list))
        end

        # @constraint(m, B_rsrv_rt <= sum(B_rt))
        # println(sum(B_feedback[bus, 1] for bus in Non_affected_Bus_list))
        # @constraint(m, B_rsrv_rt <=
        #     sum(B_rt[bus, 1] for bus in Non_affected_Bus_list))
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
                for Bus in Real_Bus_list
                    list = setdiff(Real_Bus_list, Bus)
                    @constraint(m, B_rsrv[scenario, t_real-current_time]
                        <= sum(B[scenario, bus, t_real-current_time] for bus in list))
                end

                # @constraint(m, B_rsrv[scenario, t_real-current_time]
                #      <= sum(B[scenario, :, t_real-current_time]))

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
    # println
    if sum(Pg_rt_o)<=sum(pg.mu_ct)
    ############
        Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
            +beta*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    else
        Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
            +price.lambda_ct*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    end
    # println(price.alpha_ct)
    P_cul = sum(pg.mu_ct)-sum(Pg_rt_o)
    println("curtailment solar")
    println(P_cul)
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
        pd=(pd), pg=(pg), terminate_s = (terminate_s), P_cul = (P_cul))
end

function optimal_stoach_scenario_large_solar(current_time, obj, feedback, pd, pg, price, ancillary_type, base, network, pcent_increse)
    println("===== GML - Optimization with solar issue")

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
            +(1+pcent_increse)*pg.mu[Bus,1]));
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
                    positive_scalar(icdf*sqrt(pd.sigma[Bus,t+1]+pg.sigma[Bus,t+1])
                    +(1+pcent_increse)*pg.mu[Bus,t+1]));
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

        @constraint(m, B_rsrv_rt <= sum(B_rt))
        # println(sum(B_feedback[bus, 1] for bus in Non_affected_Bus_list))
        # @constraint(m, B_rsrv_rt <=
        #     sum(B_rt[bus, 1] for bus in Non_affected_Bus_list))
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


                @constraint(m, B_rsrv[scenario, t_real-current_time]
                     <= sum(B[scenario, :, t_real-current_time]))
                # @constraint(m, B_rsrv[scenario, t_real-current_time]
                #      <= 4)
           end
       end
    end
    if ancillary_type == "10min" || ancillary_type == "30min"
        @objective(m, Min,
            fn_cost_RHC_anc_solar(delta_t,P_bus_rt,P_bus,Pg_rt,Pg,P_rsrv_rt,
            P_rsrv,price,pg,pd,beta,SN,obj,base,gen, pcent_increse))
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
    # println
    if sum(Pg_rt_o)<=sum(pg.mu_ct)
    ############
        Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
            +beta*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    else
        Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
            +price.lambda_ct*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    end
    # println(price.alpha_ct)
    P_cul = sum(pg.mu_ct)-sum(Pg_rt_o)
    println("curtailment solar")
    println(P_cul)
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
        pd=(pd), pg=(pg), terminate_s = (terminate_s), P_cul = (P_cul))
end

function optimal_stoach_scenario_large_con_solar(current_time, obj, feedback, pd, pg, price, ancillary_type, base, network, pcent_increse)
    println("===== GML - Optimization with N-1 with solar issue")

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
            +(1+pcent_increse)*pg.mu[Bus,1]));
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
                    positive_scalar(icdf*sqrt(pd.sigma[Bus,t+1]+pg.sigma[Bus,t+1])
                    +(1+pcent_increse)*pg.mu[Bus,t+1]));
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
        for Bus in Real_Bus_list
            list = setdiff(Real_Bus_list, Bus)
            # println(list)
            # println(sum(B_feedback[bus, 1] for bus in list))
            @constraint(m, B_rsrv_rt <= sum(B_rt[bus, 1] for bus in list))
        end

        # @constraint(m, B_rsrv_rt <= sum(B_rt))
        # println(sum(B_feedback[bus, 1] for bus in Non_affected_Bus_list))
        # @constraint(m, B_rsrv_rt <=
        #     sum(B_rt[bus, 1] for bus in Non_affected_Bus_list))
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
                for Bus in Real_Bus_list
                    list = setdiff(Real_Bus_list, Bus)
                    @constraint(m, B_rsrv[scenario, t_real-current_time]
                        <= sum(B[scenario, bus, t_real-current_time] for bus in list))
                end

                # @constraint(m, B_rsrv[scenario, t_real-current_time]
                #      <= sum(B[scenario, :, t_real-current_time]))

           end
       end
    end
    if ancillary_type == "10min" || ancillary_type == "30min"
        @objective(m, Min,
            fn_cost_RHC_anc_solar(delta_t,P_bus_rt,P_bus,Pg_rt,Pg,P_rsrv_rt,
            P_rsrv,price,pg,pd,beta,SN,obj,base,gen, pcent_increse))
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
    # println
    if sum(Pg_rt_o)<=sum(pg.mu_ct)
    ############
        Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
            +beta*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    else
        Cost_real = delta_t*(P_bus_rt_o[1,1]*price.lambda_ct
            +price.lambda_ct*(sum(Pg_rt_o)-sum(pg.mu_ct)))*base.MVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    end
    # println(price.alpha_ct)
    P_cul = sum(pg.mu_ct)-sum(Pg_rt_o)
    println("curtailment solar")
    println(P_cul)
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
        pd=(pd), pg=(pg), terminate_s = (terminate_s), P_cul = (P_cul))
end

function fn_cost_RHC_anc(delta_t,P_bus_rt,P_bus,Pg_rt,Pg,P_rsrv_rt,
    P_rsrv,price,pg,pd,beta,SN,obj,base,gen)

    T=obj.T;
    icdf = obj.icdf;
    gen_list = [1];
    sum_prob = sum(price.probability[1:SN])

    P_gen_sum_ct = sum(P_bus_rt[1, 1])

    Pg_ct_sum_diff = sum(Pg_rt) -
        sum(positive_array(icdf.*sqrt.(pd.sigma[:,1]+pg.sigma[:,1]))
        +pg.mu[:,1])
    Cost_Pg_ct_diff = Pg_ct_sum_diff*beta*delta_t;

    lambda_ct = sum(price.probability[scenario]/sum_prob*price.lambda_scenario[scenario,1]
    for scenario in 1:SN)

    alpha_ct = sum(price.probability[scenario]/sum_prob*price.alpha_scenario[scenario,1]
    for scenario in 1:SN)

    Cost_P_gen_sum_ct = delta_t*lambda_ct*P_gen_sum_ct;
    # println(alpha_ct)
    Revenue_P_rsrv_ct = delta_t*alpha_ct*P_rsrv_rt;

    # =======
    Cost_P_gen_sum_scenario =
        delta_t* sum(
        sum(reshape(price.probability[scenario]/sum_prob*
            price.lambda_scenario[scenario,2:T],1,T-1)*
        reshape(P_bus[scenario, 1, :],T-1,1))
        for scenario in 1:SN)

    Pg_prob_sum_scenario = sum(
        price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :])
        for scenario in 1:SN);

    Revenue_P_rsrv_scenario = sum(
        price.probability[scenario]/sum_prob*delta_t*
        reshape(P_rsrv[scenario, :], 1, T-1)*
        reshape(price.alpha_scenario[scenario, 2:end],T-1,1)
        for scenario in 1:SN);

    Pg_diff_sum_scenario = (Pg_prob_sum_scenario
        -sum(positive_array(icdf.*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])
        +pg.mu[:,:])))

    Cost_Pg_scenario_diff = Pg_diff_sum_scenario*beta*delta_t;


    println(string("    ----Case: Real-time Balancing and ", ancillary_type," Reserve Market"))

    Final_cost = ((Cost_P_gen_sum_ct+Cost_P_gen_sum_scenario
        +Cost_Pg_ct_diff+Cost_Pg_scenario_diff)*base.MVA
        -Revenue_P_rsrv_ct-Revenue_P_rsrv_scenario[1,1])
    return Final_cost
end

function fn_cost_RHC_anc_solar(delta_t,P_bus_rt,P_bus,Pg_rt,Pg,P_rsrv_rt,
    P_rsrv,price,pg,pd,beta,SN,obj,base, gen, pcent_increse)

    T=obj.T;
    icdf = obj.icdf;
    gen_list = [1];
    sum_prob = sum(price.probability[1:SN])

    P_gen_sum_ct = sum(P_bus_rt[1, 1])

    Pg_ct_sum_diff = sum(Pg_rt) -
        sum(positive_array(icdf.*sqrt.(pd.sigma[:,1]+pg.sigma[:,1]))
        +(1+ pcent_increse)*pg.mu[:,1])
    Cost_Pg_ct_diff = Pg_ct_sum_diff*beta*delta_t;

    lambda_ct = sum(price.probability[scenario]/sum_prob*price.lambda_scenario[scenario,1]
    for scenario in 1:SN)

    alpha_ct = sum(price.probability[scenario]/sum_prob*price.alpha_scenario[scenario,1]
    for scenario in 1:SN)

    Cost_P_gen_sum_ct = delta_t*lambda_ct*P_gen_sum_ct;
    # println(alpha_ct)
    Revenue_P_rsrv_ct = delta_t*alpha_ct*P_rsrv_rt;

    # =======
    Cost_P_gen_sum_scenario =
        delta_t* sum(
        sum(reshape(price.probability[scenario]/sum_prob*
            price.lambda_scenario[scenario,2:T],1,T-1)*
        reshape(P_bus[scenario, 1, :],T-1,1))
        for scenario in 1:SN)

    Pg_prob_sum_scenario = sum(
        price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :])
        for scenario in 1:SN);

    Revenue_P_rsrv_scenario = sum(
        price.probability[scenario]/sum_prob*delta_t*
        reshape(P_rsrv[scenario, :], 1, T-1)*
        reshape(price.alpha_scenario[scenario, 2:end],T-1,1)
        for scenario in 1:SN);

    Pg_diff_sum_scenario = (Pg_prob_sum_scenario
        -sum(positive_array(icdf.*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])
        +(1+pcent_increse)*pg.mu[:,:])))

    Cost_Pg_scenario_diff = Pg_diff_sum_scenario*beta*delta_t;


    println(string("    ----Case: Real-time Balancing and ", ancillary_type," Reserve Market"))
    println("with solar prediction issue")
    Final_cost = ((Cost_P_gen_sum_ct+Cost_P_gen_sum_scenario
        +Cost_Pg_ct_diff+Cost_Pg_scenario_diff)*base.MVA
        -Revenue_P_rsrv_ct-Revenue_P_rsrv_scenario[1,1])
    return Final_cost
end


function fn_cost_RHC_rt(delta_t,P_bus_rt,P_bus,Pg_rt,Pg,price,
    pg,pd,beta,SN,obj,base,gen)
    T=obj.T;
    icdf = obj.icdf;
    gen_list = [1];
    sum_prob = sum(price.probability[1:SN])



    P_gen_sum_ct = sum(P_bus_rt[gen_list[busn], 1] for busn=1:length(gen_list))
    lambda_ct = price.probability[1]/sum_prob*price.lambda_scenario[1,1]

    Pg_ct_sum_diff = sum(Pg_rt) -
        sum(positive_array(icdf.*sqrt.(pd.sigma[:,1]+pg.sigma[:,1]))
        +pg.mu[:,1])
    Cost_Pg_ct_diff = Pg_ct_sum_diff*beta*delta_t;

    if SN>1
        for scenario=2:SN
            lambda_ct =lambda_ct+price.probability[scenario]/sum_prob*price.lambda_scenario[scenario,1]
        end
    end
    Cost_P_gen_sum_ct = delta_t*lambda_ct*P_gen_sum_ct;

    # =======
    P_gen_sum_scenario = reshape(sum(P_bus[1, gen_list[busn], :] for busn=1:length(gen_list)),T-1,1)

    lambda_scenario = reshape(price.probability[1]/sum_prob*price.lambda_scenario[1,2:T],1,T-1)

    Cost_P_gen_sum_scenario = delta_t* lambda_scenario * P_gen_sum_scenario;

    Pg_prob_sum_scenario = price.probability[1]/sum_prob*sum(Pg[1, :, :])

    if SN>1
        for scenario=2:SN
            P_gen_sum_scenario = reshape(sum(P_bus[scenario, gen_list[busn], :]
                for busn=1:length(gen_list)),T-1,1)

            lambda_scenario_prob = reshape(price.probability[scenario]/
                sum_prob*price.lambda_scenario[scenario,2:T],1,T-1)

            Cost_P_gen_sum_scenario = Cost_P_gen_sum_scenario+
                delta_t* lambda_scenario_prob * P_gen_sum_scenario;

            Pg_prob_sum_scenario = Pg_prob_sum_scenario+
                price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :])
        end
    end
    Pg_diff_sum_scenario = Pg_prob_sum_scenario
        -positive_array(icdf.*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])+pg.mu[:,:])

    Cost_Pg_scenario_diff = Pg_diff_sum_scenario*beta*delta_t;


    println(string("    ----Case: Real-time Balancing"))
    return (Cost_P_gen_sum_ct+Cost_P_gen_sum_scenario[1,1]
    +Cost_Pg_ct_diff+Cost_Pg_scenario_diff)*base.MVA
end

function positive_array(Array_)
    for row=1:length(Array_[:,1])
        for column=1:length(Array_[1,:])
            if Array_[row, column]<=0
                Array_[row, column]=0;
            end
        end
    end
    return Array_
end

function positive_scalar(Scalar_)
    if Scalar_<=0
        Scalar_=0;
    end
    return Scalar_
end

function write_output_out(val_opt, P_rsrv_feed, B_5, filename)
        # write the solar file
    println("===== GML - Write Output File");
    # name=string("Time", current_time, ".csv");
    cost = val_opt.Cost_real;
    time = val_opt.time_solve;
    Pg = sum(val_opt.Pg);
    B = sum(val_opt.B);
    R = sum(val_opt.R);
    B_rsrv = sum(val_opt.B_rsrv);
    P_cul = sum(val_opt.P_cul);
    P0 = val_opt.P_bus[1,1];
    P_rsrv = P_rsrv_feed;
    status = val_opt.terminate_s;
    RT_data_feeder=hcat(cost, time, Pg, B, R, P0, P_rsrv, B_rsrv, P_cul, B_5, status)
    CSV.write(filename, DataFrame(RT_data_feeder,
        [:cost, :time, :Pg, :B, :R, :P0, :P_rsrv, :B_rsrv, :P_cul, :B_5 :status]));
    println("    ---- Finish writting files! ")
end

function write_RSRV_out(P_rsrv_feedback)
    # println(P_rsrv_feedback)
    CSV.write("P_rsrv_feedback.csv", DataFrame(reshape(P_rsrv_feedback,
        length(P_rsrv_feedback),1),
        [:P_rsrv_feedback]));
    println("    ---- Finish P_rsrv_feedback writting files! ")
end

function write_B_out(B_feedback)
    B_feedback_temp=reshape(B_feedback, length(B_feedback), 1)
    # println(B_feedback_temp)
    CSV.write("B_feedback.csv", DataFrame(B_feedback_temp,
        [:B_feedback]));
    println("    ---- Finish B_feedback writting files! ")
end

function read_RSRV_out()
    data_trace = CSV.File("P_rsrv_feedback.csv") |> DataFrame
    P_rsrv_feedback= collect(data_trace[:,Symbol("P_rsrv_feedback")])
    return reshape(P_rsrv_feedback, 1, length(P_rsrv_feedback))
end

function read_B_out()
    data_trace = CSV.File("B_feedback.csv") |> DataFrame
    B_feedback = collect(data_trace[:,Symbol("B_feedback")])
    return B_feedback
end

function read_price_data()
    filename = "../data/trace_feb28_2019.csv"
    data_trace = CSV.File(filename; dateformat="yyyy-mm-dd") |> DataFrame
    time_stamp = collect(data_trace[:,Symbol("RTD End Time Stamp")])
    ten_min_price = collect(data_trace[:,Symbol("RTD 10 Min Non Sync")])
    thirty_min_price = collect(data_trace[:,Symbol("RTD 30 Min Non Sync")])
    rt_lbmp = collect(data_trace[:, Symbol("Real Time LBMP")])
    da_lbmp = collect(data_trace[:, Symbol("Day Ahead LBMP")])
    price_raw = (timestamp = (time_stamp), RSRV_10=(ten_min_price),
        RSRV_30=(thirty_min_price), LMP_RT=(rt_lbmp), LMP_DA=(da_lbmp));
    return price_raw
end

function read_demand_data()
    filename = "../data/demand_feb_6_7.csv"
    data_trace = CSV.File(filename; dateformat="yyyy-mm-dd") |> DataFrame
    pd_da_sum= collect(data_trace[:,Symbol("da_sum")])
    # pd_feeder_2 = collect(data_trace[:,Symbol("ckt2_MW_bg276")])
    # pd_rt = hcat(pd_feeder_1, pd_feeder_2)
    pd_rt=zeros(12, 576)
    pd_da=zeros(12, 576)
    for time=1:576
        for feeder=1:12
            pd_rt[feeder, time] = data_trace[time,feeder+2];
            pd_da[feeder, time] = data_trace[time,feeder+15];
        end
    end
    pd_raw = (pd_rt = (pd_rt), pd_da_sum =(pd_da_sum), pd_da = (pd_da));
    return pd_raw
end


function read_solar_data()
    filename = "../data/Solar.csv"
    data_trace = CSV.File(filename; dateformat="yyyy-mm-dd") |> DataFrame

    # pd_feeder_2 = collect(data_trace[:,Symbol("ckt2_MW_bg276")])
    # pd_rt = hcat(pd_feeder_1, pd_feeder_2)
    pg_rt=zeros(12, 576)
    pg_da=zeros(12, 576)
    for time=1:576
        for feeder=1:12
            pg_rt[feeder, time] = data_trace[time,feeder+1];
            pg_da[feeder, time] = data_trace[(ceil(Int8, time/12))*12,feeder+1];
        end
    end
    pg_raw = (pg_rt = (pg_rt), pg_da = (pg_da));
    return pg_raw
end

function read_NY_demand_data(raw_data_mult, mult)
    norow=mult*12;
    raw_pd_one=raw_data_mult.pd_mult;
    raw_pd=repeat(raw_pd_one, mult)

    demand_unit=matread("../data/NYISO-data/normalized_demand.mat")["normalized_demand"]
    demand=zeros(norow,576)
    for row=1:norow
        demand[row,:] = reshape(
        demand_unit[:, row].*raw_pd[row], 1,576);
    end

    demand_da_unit=matread("../data/NYISO-data/normalized_demand_da.mat")["normalized_demand_da"]
    demand_da=zeros(norow,576)
    for row=1:norow
        demand_da[row,:] = reshape(
        demand_da_unit[:, row].*raw_pd[row], 1,576);
    end

    pd_raw = (pd_rt = (demand), pd_da = (demand_da));
    return pd_raw
end

function read_NY_solar_data(raw_data_mult, mult)
    norow=mult*12;
    raw_pg_one=raw_data_mult.pg_mult;
    raw_pg=repeat(raw_pg_one, mult)

    solar_unit=matread("../data/NYISO-data/normalized_solar.mat")["normalized_solar"]
    solar=zeros(norow,576)
    for row=1:norow
        solar[row,:] = reshape(
        solar_unit[:, row].*raw_pg[row], 1,576);
    end

    solar_da_unit=matread("../data/NYISO-data/normalized_solar_da.mat")["normalized_solar_da"]
    solar_da=zeros(norow,576)
    for row=1:norow
        solar_da[row,:] = reshape(
        solar_da_unit[:,row].*raw_pg[row], 1,576);
    end
    pg_raw = (pg_rt = (solar), pg_da = (solar_da));
    return pg_raw
end
