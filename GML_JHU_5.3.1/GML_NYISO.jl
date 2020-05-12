function GML_Sys_Ava_NYISO(T, pd, ancillary_type, B_cap, icdf, bus_struct,
    branch_struct, gen_struct, baseMVA)
    println("===== GML - Boundaries Buildup");
    ###############################################################################
    T=288;
    NoBus = length(bus_struct.baseKV)
    NoBr = length(branch_struct.fbus)
    Qf_max = zeros(NoBus,T)
    Qf_min = zeros(NoBus,T)
    for bus = 1:NoBus
        # println(bus)
        if bus_struct.type[bus] == 1
            Qf_max[bus, :]=0.05*positive_array(pd.traj[bus, :]);
            if sum(Qf_max[bus,:])==0
                Qf_max[bus, :] = 0.1*ones(1,T)
            end
            # println(sum(Qf_max[bus,:]))
            Qf_min = -Qf_max;
            # println(sum(Qf_min[bus,:]))
        else
            gen_id = findall(gen->gen==bus, gen_struct.bus)[1];
            # println(gen_id)
            Qf_max[bus, :]=ones(1,T)*gen_struct.Qmax[gen_id];
            # println(sum(Qf_max[bus,:]))
            Qf_min[bus, :]=ones(1,T)*gen_struct.Qmin[gen_id]
            # println(sum(Qf_min[bus,:]))
        end
    end


    # minimum solar
    Pg_min = zeros(NoBus, T);
    B_max = reshape(bus_struct.frac*B_cap, NoBus, 1)*ones(1, T);
    # println(sum(B_max)/288)
    B_min = zeros(NoBus, T);
    R_rate = 1/3;
    R_max = R_rate*B_max/baseMVA;
    R_min = -R_max;
    # println(R_max[:,1])
    # println(R_min[:,1])
    W = zeros(NoBus, T);
    #
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
    #
    V_min = 0.94;
    V_max = 1.06;
    r=branch_struct.r
    # println(r)
    x=branch_struct.x
    # println(x)
    # # # penalty
    beta=-15;

    C_ind = zeros(NoBus, NoBr)
    for branch = 1:NoBr
        C_ind[Int(branch_struct.fbus[branch]), branch]=1;
        C_ind[Int(branch_struct.tbus[branch]), branch]=-1;
    end
    # println(C_ind)
    println("===== GML - finish building boundaries =====")
    obj=(T=(T), icdf = (icdf), Qf_max=(Qf_max), Qf_min=(Qf_min), Pg_min=(Pg_min),
    B_max=(B_max), B_min=(B_min), R_max=(R_max), R_min=(R_min),
    tau=(tau), P_rsrv_min=(P_rsrv_min), k=(k),V_max=(V_max), V_min=(V_min),
    r=(r), x=(x), beta=(beta), C_ind=(C_ind));

    return obj
end

function optimal_NYISO(SN, t, obj, ancillary_type, baseMVA,
    feedback, pd, pg, price, bus_struct, branch_struct, gen_struct);
    Q_gamma=0.01;
    println("===== GML - Optimization ")
    T=obj.T;
    Pg_min = obj.Pg_min; # minimum power generation

    Qf_max = obj.Qf_max; # maximum reactive power generation
    Qf_min = obj.Qf_min; # minimum eactive power generation

    B_max = obj.B_max; # maximum storage level
    B_min = obj.B_min; # minimum storage level

    R_max = obj.R_max; # the maximum discharge rate
    R_min = obj.R_min; # the minimum discharge rate

    delta_t = 1/12;; # time interval
    icdf = obj.icdf;
    C_ind = obj.C_ind;

    # P_rsrv_max = obj.P_rsrv_max; # The maximum ancillary power

    P_rsrv_min = obj.P_rsrv_min; # The minimum ancillary power

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

    NoBus = length(bus_struct.baseKV)
    NoBr = length(branch_struct.fbus)
    NoGen = length(gen_struct.Pmax)
    # NoGen=length(gen_data.id)

    m = Model(with_optimizer(Mosek.Optimizer, QUIET=true,
    MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-3,
    MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-3))
    # define the real-time variables

    # Bus level
    @variable(m, Pg_rt[1:NoBus, 1])
    # @variable(m, Qf_rt[1:NoBus, 1])
    @variable(m, B_rt[1:NoBus, 1])
    @variable(m, R_rt[1:NoBus, 1])

    @variable(m, P_bus_rt[1:NoBus, 1])
    @variable(m, Q_bus_rt[1:NoBus, 1])
    @variable(m, v_rt[1:NoBus, 1])
    # Branch level
    @variable(m, P_br_rt[1:NoBr, 1])
    @variable(m, Q_br_rt[1:NoBr, 1])

    @variable(m, P_gen_rt[1:NoGen, 1])
    @variable(m, Q_gen_rt[1:NoGen, 1])


    # println(" ---- Real Time Constraint Buildup ")
    for bus=1:NoBus
        # box constraints on Solar (Pg), Battery discharge (R)
        @constraint(m, 0<=Pg_rt[bus,1]);
        @constraint(m, Pg_rt[bus,1]<=
            positive_scalar(icdf*sqrt(pd.sigma[bus,1]+pg.sigma[bus,1])
            +pg.mu[bus,1]));
        @constraint(m, R_min[bus,1]<= R_rt[bus,1]);
        @constraint(m, R_rt[bus,1]<= R_max[bus,1]);
        # Initial battery SOC
        @constraint(m, B_rt[bus,1]==B_feedback[bus,1]);

        # SOC constrains on real and reactive power on bus
        maxS_rt = sqrt(1+Q_gamma^2)*(Pd[bus,1]+R_max[bus,1]-
        positive_scalar(
        icdf*sqrt(pd.sigma[bus,1]+pg.sigma[bus,1])));
        if maxS_rt>=0
            @constraint(m,
            [maxS_rt,
            P_bus_rt[bus,1], Q_bus_rt[bus,1]] in SecondOrderCone())
        else
            @constraint(m,
            [-maxS_rt,
            P_bus_rt[bus,1], Q_bus_rt[bus,1]] in SecondOrderCone())
        end
        # a list of branchs that FROM the bus of interest
        sub_branch_list = findall(one->one==1, C_ind[bus,:]) # Out
        # a list of branchs that POINT TO the bus of interest
        add_branch_list = findall(minusone->minusone==-1, C_ind[bus,:]) #In

        # real power on bus is demand minus solar and discharge
        # AKA power injection
        @constraint(m, P_bus_rt[bus,1]==
        Pd[bus,1]
        -(Pg_rt[bus,1]+R_rt[bus,1]));

        if bus_struct.type[bus]==1 # non-generator bus
            # box constraint on voltage
            @constraint(m, V_min^2<= v_rt[bus,1]);
            @constraint(m, v_rt[bus,1]<= V_max^2);

            # Power Balance Equations
            # Power injection = Power Flow In - Power Flow Out
            if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    sum(P_br_rt[branch,1] for branch in add_branch_list)
                    -sum(P_br_rt[branch,1] for branch in sub_branch_list)
                    );
                @constraint(m, Q_bus_rt[bus,1]==
                    sum(Q_br_rt[branch,1] for branch in add_branch_list)
                    -sum(Q_br_rt[branch,1] for branch in sub_branch_list));
            elseif ~isempty(add_branch_list) && isempty(sub_branch_list)

                @constraint(m, P_bus_rt[bus,1]==
                    sum(P_br_rt[branch,1] for branch in add_branch_list)
                    );
                @constraint(m, Q_bus_rt[bus,1]==
                    sum(Q_br_rt[branch,1] for branch in add_branch_list));
            elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    -sum(P_br_rt[branch,1] for branch in sub_branch_list)
                    );
                @constraint(m, Q_bus_rt[bus,1]==
                    -sum(Q_br_rt[branch,1] for branch in sub_branch_list));
            end
        else # generator bus
            ########
            # @constraint(m, v_rt[bus,1]== 1);
            # only one generator, and is regraded as the slack bus
            ########
            @constraint(m, V_min^2<= v_rt[bus,1]);
            @constraint(m, v_rt[bus,1]<= V_max^2);

            # identify the id of the generator that connects to the bus of interest
            gen_id = findall(bus_id ->bus_id == bus, gen_struct.bus);
            # Power Balance Equations
            # Power injection = Power Flow In - Power Flow Out + Power Generatered
            if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    sum(P_br_rt[branch,1] for branch in add_branch_list)
                    -sum(P_br_rt[branch,1] for branch in sub_branch_list)
                    +P_gen_rt[gen_id[1], 1]);
                @constraint(m, Q_bus_rt[bus,1]==
                    sum(Q_br_rt[branch,1] for branch in add_branch_list)
                    -sum(Q_br_rt[branch,1] for branch in sub_branch_list)
                    +Q_gen_rt[gen_id[1], 1]);
            elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    sum(P_br_rt[branch,1] for branch in add_branch_list)
                    +P_gen_rt[gen_id[1], 1]);
                @constraint(m, Q_bus_rt[bus,1]==
                    sum(Q_br_rt[branch,1] for branch in add_branch_list)
                    +Q_gen_rt[gen_id[1], 1]);
            elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    -sum(P_br_rt[branch,1] for branch in sub_branch_list)
                    +P_gen_rt[gen_id[1], 1]);
                @constraint(m, Q_bus_rt[bus,1]==
                    -sum(Q_br_rt[branch,1] for branch in sub_branch_list)
                    +Q_gen_rt[gen_id[1], 1]);
            end
        end
    end
    # voltage constraint for LinDistFlow
    for branch =1:NoBr
        fbus = branch_struct.fbus[branch];
        tbus = branch_struct.tbus[branch];
        @constraint(m, v_rt[fbus]-v_rt[tbus]==
            2*(r[branch]*P_br_rt[branch]+x[branch]*Q_br_rt[branch]))
    end
    # box constraint on generator
    # NoGen
    for Gen=1:NoGen
        @constraint(m,
            P_gen_rt[Gen,1]<=gen_struct.Pmax[Gen]/baseMVA)
        @constraint(m,
            P_gen_rt[Gen,1]>=gen_struct.Pmin[Gen]/baseMVA)
        @constraint(m,
            Q_gen_rt[Gen,1]<=gen_struct.Qmax[Gen]/baseMVA)
        @constraint(m,
            Q_gen_rt[Gen,1]>=gen_struct.Qmin[Gen]/baseMVA)
    end

    # bus
    @variable(m, Pg[1:SN, 1:NoBus, 1:T-1]); # the real power output
    @variable(m, B[1:SN, 1:NoBus, 1:T-1]); # the storage
    @variable(m, R[1:SN, 1:NoBus, 1:T-1]);# the charge/discharge rate
    @variable(m, P_bus[1:SN, 1:NoBus, 1:T-1])
    @variable(m, Q_bus[1:SN, 1:NoBus, 1:T-1])
    @variable(m, v[1:SN, 1:NoBus, 1:T-1])

    # Branch level
    @variable(m, P_br[1:SN,1:NoBr, 1:T-1])
    @variable(m, Q_br[1:SN,1:NoBr, 1:T-1])
    # Gen Level
    @variable(m, P_gen[1:SN, 1:NoGen, 1:T-1])
    @variable(m, Q_gen[1:SN, 1:NoGen, 1:T-1])

    # for different scenario
    for scenario = 1:SN
        # for different time at prediction horizion
        for t=1:T-1
            # for different bus
            for bus=1:NoBus
                # box constraints on Solar (Pg), Battery discharge (R)
                @constraint(m, Pg_min[bus, t]<=Pg[scenario, bus ,t]);
                @constraint(m, Pg[scenario, bus,t]<=
                    positive_scalar(
                    icdf*sqrt(pd.sigma[bus,t+1]+pg.sigma[bus,t+1])+pg.mu[bus,t+1])
                    );
                @constraint(m, R_min[bus,t+1] <= R[scenario,bus,t]);
                @constraint(m, R[scenario,bus,t]<= R_max[bus,t+1]);
                # battery box constraint on SOC
                @constraint(m, B_min[bus,t+1] <= B[scenario,bus,t]);
                @constraint(m, B[scenario,bus,t]<= B_max[bus,t+1]);
                # battery SOC dynamics
                if t==1
                    @constraint(m, B[scenario,bus,1] == B_rt[bus,1]
                        -delta_t*(R_rt[bus,1]*baseMVA))
                else
                    @constraint(m, B[scenario,bus,t] ==
                        B[scenario,bus,t-1] - R[scenario,bus,t-1]*baseMVA*delta_t)
                end
                # real power on bus is demand minus solar and discharge
                # AKA power injection
                @constraint(m,  P_bus[scenario,bus,t]==
                Pd[bus,t+1]
                -Pg[scenario,bus,t]-R[scenario, bus,t]);

                # SOC constrains on real and reactive power on bus
                maxS = sqrt(1+Q_gamma^2)*(Pd[bus,t+1]-positive_scalar(
                icdf*sqrt(pd.sigma[bus,t+1]+pg.sigma[bus,t+1])+pg.mu[bus,t+1])
                -R_max[bus,t+1]);
                if maxS >=0
                    @constraint(m,
                    [maxS, P_bus[scenario,bus,t], Q_bus[scenario,bus,t]]
                     in SecondOrderCone())
                else
                    @constraint(m,
                    [-maxS, P_bus[scenario,bus,t], Q_bus[scenario,bus,t]]
                     in SecondOrderCone())
                end

                # a list of branchs that FROM the bus of interest
                sub_branch_list = findall(one->one==1, C_ind[bus,:]) # out
                # a list of branchs that POINT TO the bus of interest
                add_branch_list = findall(minusone->minusone==-1, C_ind[bus,:]) # In

                if bus_struct.type[bus]==1 # non-generator bus
                    # box constraint on bus
                    @constraint(m, V_min^2<= v[scenario,bus,t]);
                    @constraint(m, v[scenario,bus,t]<= V_max^2);
                    # Power Balance Equations
                    # Power injection = Power Flow In - Power Flow Out
                    if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            sum(P_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
                            );
                        @constraint(m, Q_bus[scenario,bus,t]==
                            sum(Q_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
                            );
                    elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            sum(P_br[scenario,branch,t] for branch in add_branch_list)
                            );
                        @constraint(m, Q_bus[scenario,bus,t]==
                            sum(Q_br[scenario,branch,t] for branch in add_branch_list)
                            );
                    elseif isempty(add_branch_list) && ~isempty(sub_branch_list)

                        @constraint(m, P_bus[scenario,bus,t]==
                            -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
                            );
                        @constraint(m, Q_bus[scenario,bus,t]==
                            -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
                            );
                    end
                else # genertor bus

                    @constraint(m, V_min^2<= v[scenario,bus,t]);
                    @constraint(m, v[scenario,bus,t]<= V_max^2);
                    gen_id = findall(bus_id ->bus_id == bus, gen_struct.bus);
                    # @constraint(m, v[scenario,bus,t]==1);

                    if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            sum(P_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
                            +P_gen[scenario,gen_id[1],t]);
                        @constraint(m, Q_bus[scenario,bus,t]==
                            sum(Q_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
                            +Q_gen[scenario,gen_id[1],t]);
                    elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            sum(P_br[scenario,branch,t] for branch in add_branch_list)
                            +P_gen[scenario,gen_id[1],t]);
                        @constraint(m, Q_bus[scenario,bus,t]==
                            sum(Q_br[scenario,branch,t] for branch in add_branch_list)
                            +Q_gen[scenario,gen_id[1],t]);
                    elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
                            +P_gen[scenario,gen_id[1],t]);
                        @constraint(m, Q_bus[scenario,bus,t]==
                            -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
                            +Q_gen[scenario,gen_id[1],t]);
                    end
                end
            end

            for branch =1:NoBr
                fbus = Int(branch_struct.fbus[branch]);
                tbus = Int(branch_struct.tbus[branch]);
                @constraint(m, v[scenario,fbus,t]-v[scenario,tbus,t]==
                    2*(r[branch]*P_br[scenario,branch,t]
                    +x[branch]*Q_br[scenario,branch,t]))
            end

            for Gen=1:NoGen
                @constraint(m,
                    P_gen[scenario,Gen,t]<=gen_struct.Pmax[Gen]/baseMVA)
                @constraint(m,
                    P_gen[scenario,Gen,t]>=gen_struct.Pmin[Gen]/baseMVA)
                @constraint(m,
                    Q_gen[scenario,Gen,t]<=gen_struct.Qmax[Gen]/baseMVA)
                @constraint(m,
                    Q_gen[scenario,Gen,t]>=gen_struct.Qmin[Gen]/baseMVA)
            end
        end
    end


    if ancillary_type == "10min" || ancillary_type == "30min"
            # println(Real_Bus_list)
            # RT
            @variable(m, P_rsrv_rt)
            @variable(m, B_rsrv_rt)

            @constraint(m, P_rsrv_rt>=0)
            @constraint(m, B_rsrv_rt>=0)

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
            # for Bus in Real_Bus_list
            #     list = setdiff(Real_Bus_list, Bus)
            #     # println(list)
            #     # println(sum(B_feedback[bus, 1] for bus in list))
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
                    # for Bus in Real_Bus_list
                    #     list = setdiff(Real_Bus_list, Bus)
                    #     @constraint(m, B_rsrv[scenario, t_real-current_time]
                    #         <= sum(B[scenario, bus, t_real-current_time] for bus in list))
                    # end

                    @constraint(m, B_rsrv[scenario, t_real-current_time]
                         <= sum(B[scenario, :, t_real-current_time]))

               end
           end
        end
    # #

    if ancillary_type == "10min" || ancillary_type == "30min"
        @objective(m, Min,
            fn_cost_RHC_anc(delta_t, P_gen, P_gen_rt, Pg_rt,Pg, P_rsrv_rt,
            P_rsrv,price,pg,pd,beta,SN,obj, gen_struct, bus_struct, branch_struct))
    else
        @objective(m, Min,
            fn_cost_RHC_rt(delta_t, P_gen, P_gen_rt, Pg_rt,Pg, price,
                pg, pd, beta, SN, obj, gen_struct, bus_struct, branch_struct, baseMVA))
    end
    #
    status=optimize!(m);
    println(string("    ----", termination_status(m)))
    terminate_s = termination_status(m);
    # println(MOI.PrimalStatus())
    # println(MOI.DualStatus())
    cost_o = JuMP.objective_value(m);
    time_solve=MOI.get(m, MOI.SolveTime());
    println(string("    ----Solve Time: ", time_solve))
    println(string("    ----Optimal Cost of the RHC: ", cost_o))

    println(string("    ----aggregate demand: ", baseMVA*sum(Pd[:,1])))
    ###############
    R_rt_o=JuMP.value.(R_rt)
    R_o=JuMP.value.(R)
    R_sum_ct = sum(R_rt_o[:, 1]*baseMVA)
    println(string("    ----aggregate battery disharge: ", R_sum_ct))
    R_traj = hcat(R_rt_o, R_o[1,:,:])
    ###############
    P_br_rt_o=JuMP.value.(P_br_rt)
    P_br_o=JuMP.value.(P_br)
    P_br_traj = hcat(P_br_rt_o, P_br_o[1,:,:])

    ###############
    Q_br_rt_o=JuMP.value.(Q_br_rt)
    Q_br_o=JuMP.value.(Q_br)
    Q_br_traj = hcat(Q_br_rt_o, Q_br_o[1,:,:])
    ###############
    v_rt_o=JuMP.value.(v_rt)
    v_o=JuMP.value.(v)
    v_traj = hcat(v_rt_o, v_o[1,:,:])
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
    Q_bus_rt_o=JuMP.value.(Q_bus_rt)
    Q_bus_o=JuMP.value.(Q_bus)
    Q_bus_traj = hcat(Q_bus_rt_o, Q_bus_o[1,:,:])
    println(string("    ----aggregate reactive on bus: ",
        sum(Q_bus_traj)))

    ############
    P_gen_rt_o=JuMP.value.(P_gen_rt)
    P_gen_o=JuMP.value.(P_gen)
    P_gen_traj = hcat(P_gen_rt_o, P_gen_o[1,:,:])
    P_gen_sum_traj = hcat(sum(P_gen_rt_o), sum(P_gen_o[1,:,:], dims=1))

    P_gen_sum_ct = sum(P_gen_rt_o[gen, 1]*baseMVA for gen in 1:NoGen)
    println(string("    ----aggregate generation: ", P_gen_sum_ct))


    ############
    Q_gen_rt_o = JuMP.value.(Q_gen_rt)
    Q_gen_o = JuMP.value.(Q_gen)
    Q_gen_traj = hcat(Q_gen_rt_o, Q_gen_o[1,:,:])
    Q_gen_sum_traj = hcat(sum(Q_gen_rt_o), sum(Q_gen_o[1,:,:], dims=1))
    Q_gen_sum_ct = sum(Q_gen_rt_o[gen, 1]*baseMVA for gen in 1:NoGen)
    println(string("    ----aggregate reactive generation: ",
        sum(Q_gen_sum_ct)))
    ###########
    # println("print voltage constraint")
    # for branch =15:16
    #     fbus = branch_struct.fbus[branch];
    #     tbus = branch_struct.tbus[branch];
    #     println(v_rt_o[fbus]-v_rt_o[tbus])
    #     # println(v_rt_o[fbus]-v_rt_o[tbus]-
    #     #     2*(r[branch]*P_br_rt_o[branch,1]))
    # end
    # for t=1:287
    #     for branch =15:16
    #         fbus = Int(branch_struct.fbus[branch]);
    #         tbus = Int(branch_struct.tbus[branch]);
    #         println(v_o[1,fbus,t]-v_o[1,tbus,t])
    #         # println(v_o[1,fbus,t]-v_o[1,tbus,t]-
    #         #     2*(r[branch]*P_br_o[1,branch,t]
    #         #     ))
    #     end
    # end

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
        B_rsrv_total = zeros(1,T);
    end
    # println
    if sum(Pg_rt_o)<=sum(pg.mu_ct)
    ############
        Cost_real = delta_t*(sum(P_gen_rt_o)*price.lambda_ct
            +beta*(sum(Pg_rt_o)-sum(pg.mu_ct)))*baseMVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    else
        Cost_real = delta_t*(sum(P_gen_rt_o)*price.lambda_ct
            +price.lambda_ct*(sum(Pg_rt_o)-sum(pg.mu_ct)))*baseMVA-
            delta_t*P_rsrv_rt_o*price.alpha_ct;
    end
   # println(price.alpha_ct)
    P_cul = sum(pg.mu_ct)-sum(Pg_rt_o)
    println(string("    ----curtailment solar: ", P_cul))
    alpha_1 = price.alpha_scenario[1,:]
    lambda_1 = price.lambda_scenario[1,:]

    println(string("    ----Optimzal cost at this instance: ", Cost_real))
    pg_upper = hcat(
    icdf*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])+pg.mu[:,:])

    val_opt = (R=(R_rt_o), B=(B_rt_o), Pg=(Pg_rt_o), P_gen=(P_gen_rt_o), Q_gen=(Q_gen_rt_o),
        P_rsrv=(P_rsrv_rt_o), B_rsrv=(B_rsrv_rt_o), P_bus=(P_bus_rt_o), Q_bus=(Q_bus_rt_o),
        v=(v_rt_o), P_br=(P_br_rt_o), Q_br=(Q_br_rt_o),
        Cost_real=(Cost_real), time_solve=(time_solve),
        P_rsrv_total=(P_rsrv_total), B_rsrv_total=(B_rsrv_total), alpha_1=(alpha_1),
        R_traj = (R_traj), B_traj=(B_traj), Pg_traj=(Pg_traj), P_gen_traj=(P_gen_traj),
        Q_gen_traj=(Q_gen_traj),
        P0_traj = (P_gen_traj), P_bus_traj=(P_bus_traj), Q_bus_traj=(Q_bus_traj),
        pg_upper=(pg_upper),
        lambda_1=(lambda_1), pd=(pd), pg=(pg), v_traj=(v_traj), P_br_traj=(P_br_traj),
        Q_br_traj=(Q_br_traj), terminate_s = (terminate_s), P_cul = (P_cul))


    return val_opt
end


function fn_cost_RHC_rt(delta_t, P_gen, P_gen_rt, Pg_rt,Pg, price,
    pg, pd, beta, SN, obj, gen_struct, bus_struct, branch_struct, baseMVA)

    T=obj.T;
    icdf = obj.icdf;


    sum_prob = sum(price.probability[1:SN]);
    NoGen = length(gen_struct.Pmax);

    ################### Generation Cost ###################
    # At current time
    lambda_ct = sum(price.probability[scenario]/sum_prob*price.lambda_scenario[scenario,1]
    for scenario in 1:SN)
    P_gen_sum_ct = sum((P_gen_rt[gen, 1])*baseMVA for gen in 1:NoGen);
    Cost_P_gen_sum_ct = delta_t*lambda_ct*P_gen_sum_ct;
    # At future time
    Cost_P_gen_sum_scenario =
        delta_t*sum(
        sum(reshape(price.probability[scenario]/sum_prob*
            price.lambda_scenario[scenario,2:T],1,T-1)*
        reshape(
        sum(P_gen[scenario, gen, :] for gen in 1:NoGen),
        T-1,1))
        for scenario in 1:SN)*baseMVA;

    ################### FOL Cost (Solar penalty)###################

    Pg_ct_sum_diff = sum(Pg_rt) -
        sum(positive_array(icdf.*sqrt.(pd.sigma[:,1]+pg.sigma[:,1]))
        +pg.mu[:,1]);
    Cost_Pg_ct_diff = Pg_ct_sum_diff*beta*delta_t*baseMVA;

    Pg_prob_sum_scenario = sum(
        price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :])
        for scenario in 1:SN);
    Pg_diff_sum_scenario = (Pg_prob_sum_scenario
        -sum(positive_array(icdf.*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])+pg.mu[:,:])));
    Cost_Pg_scenario_diff = Pg_diff_sum_scenario*beta*delta_t*baseMVA;


    Final_cost = (Cost_P_gen_sum_ct+Cost_P_gen_sum_scenario
        +Cost_Pg_ct_diff+Cost_Pg_scenario_diff)

    return Final_cost
end

function fn_cost_RHC_anc(delta_t, P_gen, P_gen_rt, Pg_rt,Pg, P_rsrv_rt,
    P_rsrv,price,pg,pd,beta,SN,obj, gen_struct, bus_struct, branch_struct)

    T=obj.T;
    icdf = obj.icdf;


    sum_prob = sum(price.probability[1:SN]);
    NoGen = length(gen_struct.Pmax);

    ################### Generation Cost ###################
    # At current time
    lambda_ct = sum(price.probability[scenario]/sum_prob*price.lambda_scenario[scenario,1]
    for scenario in 1:SN)
    P_gen_sum_ct = sum((P_gen_rt[gen, 1])*baseMVA for gen in 1:NoGen);
    Cost_P_gen_sum_ct = delta_t*lambda_ct*P_gen_sum_ct;
    # At future time
    Cost_P_gen_sum_scenario =
        delta_t*sum(
        sum(reshape(price.probability[scenario]/sum_prob*
            price.lambda_scenario[scenario,2:T],1,T-1)*
        reshape(
        sum(P_gen[scenario, gen, :] for gen in 1:NoGen),
        T-1,1))
        for scenario in 1:SN)*baseMVA;

    ################### FOL Cost (Solar penalty)###################

    Pg_ct_sum_diff = sum(Pg_rt) -
        sum(positive_array(icdf.*sqrt.(pd.sigma[:,1]+pg.sigma[:,1]))
        +pg.mu[:,1]);
    Cost_Pg_ct_diff = Pg_ct_sum_diff*beta*delta_t*baseMVA;

    Pg_prob_sum_scenario = sum(
        price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :])
        for scenario in 1:SN);
    Pg_diff_sum_scenario = (Pg_prob_sum_scenario
        -sum(positive_array(icdf.*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])+pg.mu[:,:])));
    Cost_Pg_scenario_diff = Pg_diff_sum_scenario*beta*delta_t*baseMVA;

    ################### Revenue from ancillary ###################
    alpha_ct = sum(price.probability[scenario]/sum_prob*price.alpha_scenario[scenario,1]
        for scenario in 1:SN);
    Revenue_P_rsrv_ct = delta_t*alpha_ct*P_rsrv_rt;
    Revenue_P_rsrv_scenario = sum(
        sum(
        price.probability[scenario]/sum_prob*delta_t*
        reshape(P_rsrv[scenario, :], 1, T-1)*
        reshape(price.alpha_scenario[scenario, 2:end],T-1,1)
        for scenario in 1:SN)
        );

    Final_cost = (Cost_P_gen_sum_ct+Cost_P_gen_sum_scenario
        +Cost_Pg_ct_diff+Cost_Pg_scenario_diff
        -Revenue_P_rsrv_ct-Revenue_P_rsrv_scenario)
    # Final_cost = (Cost_P_gen_sum_ct)
    return Final_cost
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

###############################
function read_NY_demand_data(bus)
    # P=shunt.P;
    NoBus = length(bus.baseKV)
    demand_unit=matread("../data/NYISO-data/normalized_demand.mat")["normalized_demand"]
    n = size(demand_unit)
    demand_bus = zeros(n[2], n[1])
    for i = 1 : n[1]
        demand_bus[:,i] = reshape(demand_unit[i,:],(n[2], 1))
    end
    frac = bus.frac

    for i = 1 : NoBus
        demand_bus[i,:] = demand_bus[i,:] .* frac[i]
    end

    demand_bus = demand_bus[1:NoBus,:]

    demand_da_bus = zeros(NoBus, n[1])

    # demand_da_bus

    # demand=zeros(length(P),576)
    # for n_shunt=1:length(P)
    #     demand[n_shunt,:] = reshape(
    #     demand_unit[:, n_shunt].*(576).*P[n_shunt].*0.7, 1,576);
    # end
    #
    # demand_da_unit=matread("../data/NYISO-data/normalized_demand_da.mat")["normalized_demand_da"]
    # demand_da=zeros(length(P),576)
    # for n_shunt=1:length(P)
    #     demand_da[n_shunt,:] = reshape(
    #     demand_da_unit[:, n_shunt].*(576).*P[n_shunt].*0.7, 1,576);
    # end
    #
    # demand_bus = zeros(NoBus, 576)
    # demand_da_bus = zeros(NoBus, 576)
    #
    # for bus = 1:NoBus
    #     shunt_list = findall(inbus->inbus==bus, shunt.find_bus)
    #     if isempty(shunt_list)
    #         demand_bus[bus, :]=zeros(1, 576)
    #         demand_da_bus[bus, :]=zeros(1, 576)
    #     else
    #         demand_bus[bus, :] = sum(demand[shunt, :] for shunt in shunt_list)
    #         demand_da_bus[bus, :] = sum(demand_da[shunt, :] for shunt in shunt_list)
    #     end
    # end
    #
    pd_raw = (pd_rt = (demand_bus), pd_da = (demand_da_bus));
    return pd_raw
end

###################################
function read_NY_solar_data(bus)
    # P=shunt.P;
    NoBus = length(bus.baseKV);
    solar_unit=matread("../data/NYISO-data/normalized_solar.mat")["normalized_solar"]
    n = size(solar_unit)
    solar_bus = zeros(n[2], n[1])
    for i = 1 : n[1]
        solar_bus[:,i] = reshape(solar_unit[i,:],(n[2], 1))
    end
    frac = bus.frac

    for i = 1 : NoBus
        solar_bus[i,:] = solar_bus[i,:] .* frac[i]
    end

    solar_bus = solar_bus[1:NoBus,:]

    solar_da_bus = zeros(NoBus, n[1])

    # solar=zeros(length(P),576)
    # for n_shunt=1:length(P)
    #     solar[n_shunt,:] = reshape(
    #     solar_unit[:, n_shunt].*(576).*P[n_shunt].*0.7, 1,576);
    # end
    #
    # solar_da_unit=matread("../data/NYISO-data/normalized_solar_da.mat")["normalized_solar_da"]
    # solar_da=zeros(length(P),576)
    # for n_shunt=1:length(P)
    #     solar_da[n_shunt,:] = reshape(
    #     solar_da_unit[:, n_shunt].*(576).*P[n_shunt].*0.7, 1,576);
    # end
    # solar_bus = zeros(NoBus, 576)
    # solar_da_bus = zeros(NoBus, 576)
    # for bus = 1:NoBus
    #     shunt_list = findall(inbus->inbus==bus, shunt.find_bus)
    #     if isempty(shunt_list)
    #         solar_bus[bus, :]=zeros(1, 576)
    #         solar_da_bus[bus, :]=zeros(1, 576)
    #     else
    #         solar_bus[bus, :] = sum(solar[shunt, :] for shunt in shunt_list)
    #         solar_da_bus[bus, :] = sum(solar_da[shunt, :] for shunt in shunt_list)
    #     end
    # end
    pg_raw = (pg_rt = (solar_bus), pg_da = (solar_da_bus));
    return pg_raw
end


function write_branch_real_output(val_opt)
    output = hcat(val_opt.P_br_traj)
    CSV.write("branch_flow_real_power_output.csv", DataFrame(output));
    println("    ---- Finish branch flow real files! ")
end

function write_branch_reactive_output(val_opt)
    output = hcat(val_opt.Q_br_traj)
    CSV.write("branch_flow_reactive_power_output.csv", DataFrame(output));
    println("    ---- Finish branch flow reactive files! ")
end

function write_bus_real_output(val_opt)
    output = hcat(val_opt.P_bus_traj)
    CSV.write("bus_real_output.csv", DataFrame(output));
    println("    ---- Finish bus real files! ")
end

function write_bus_reactive_output(val_opt)
    output = hcat(val_opt.Q_bus_traj)
    CSV.write("bus_reactive_output.csv", DataFrame(output));
    println("    ---- Finish bus reactive files! ")
end

function write_bus_voltage_output(val_opt)
    output = hcat(sqrt.(val_opt.v_traj))
    CSV.write("bus_voltage_output.csv", DataFrame(output));
    println("    ---- Finish bus voltage files! ")
end

function write_generator_real_output(val_opt)
    output = hcat(val_opt.P_gen_traj)
    CSV.write("gen_real_output.csv", DataFrame(output));
    println("    ---- Finish bus files! ")
end

function write_generator_reactive_output(val_opt)
    output = hcat(val_opt.Q_gen_traj)
    CSV.write("gen_reactive_output.csv", DataFrame(output));
    println("    ---- Finish bus files! ")
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

# function optimal_NYISO_SOCP(SN, t, obj, ancillary_type, baseMVA,
#     feedback, pd, pg, price, bus_struct, branch_struct, gen_struct);
#     Q_gamma=0.05;
#     println("===== GML - Optimization ")
#     T=obj.T;
#     Pg_min = obj.Pg_min; # minimum power generation
#
#     Qf_max = obj.Qf_max; # maximum reactive power generation
#     Qf_min = obj.Qf_min; # minimum eactive power generation
#
#     B_max = obj.B_max; # maximum storage level
#     B_min = obj.B_min; # minimum storage level
#
#     R_max = obj.R_max; # the maximum discharge rate
#     R_min = obj.R_min; # the minimum discharge rate
#
#     delta_t = 1/12;; # time interval
#     icdf = obj.icdf;
#     C_ind = obj.C_ind;
#
#     # P_rsrv_max = obj.P_rsrv_max; # The maximum ancillary power
#
#     P_rsrv_min = obj.P_rsrv_min; # The minimum ancillary power
#
#     V_min=obj.V_min;
#     V_max=obj.V_max;
#
#     beta=obj.beta; # price in cost function
#     tau = obj.tau;
#
#     r=obj.r; # the resistance
#     x=obj.x; # the reactance
#     k=obj.k; # time with ancellary
#
#     Pd = pd.traj;
#
#     B_feedback = feedback.B_feedback;
#     P_rsrv_feedback = feedback.P_rsrv_feedback;
#     # println(size(feedback.P_rsrv_feedback))
#
#     NoBus = length(bus_struct.baseKV)
#     NoBr = length(branch_struct.fbus)
#     NoGen = length(gen_struct.Pmax)
#     # NoGen=length(gen_data.id)
#
#     m = Model(with_optimizer(Mosek.Optimizer, QUIET=true,
#     MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-3,
#     MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-3))
#     # define the real-time variables
#
#     # Bus level
#     @variable(m, Pg_rt[1:NoBus, 1])
#     # @variable(m, Qf_rt[1:NoBus, 1])
#     @variable(m, B_rt[1:NoBus, 1])
#     @variable(m, R_rt[1:NoBus, 1])
#
#     @variable(m, P_bus_rt[1:NoBus, 1])
#     @variable(m, Q_bus_rt[1:NoBus, 1])
#     @variable(m, v_rt[1:NoBus, 1])
#
#     # Branch level
#     @variable(m, P_br_rt[1:NoBr, 1])
#     @variable(m, Q_br_rt[1:NoBr, 1])
#     @variable(m, l_br_rt[1:NoBr, 1])
#
#     @variable(m, P_gen_rt[1:NoGen, 1])
#     @variable(m, Q_gen_rt[1:NoGen, 1])
#
#
#     # println(" ---- Real Time Constraint Buildup ")
#     for bus=1:NoBus
#         # box constraints on Solar (Pg), Battery discharge (R)
#         @constraint(m, 0<=Pg_rt[bus,1]);
#         @constraint(m, Pg_rt[bus,1]<=
#             positive_scalar(icdf*sqrt(pd.sigma[bus,1]+pg.sigma[bus,1])
#             +pg.mu[bus,1]));
#         @constraint(m, R_min[bus,1]<= R_rt[bus,1]);
#         @constraint(m, R_rt[bus,1]<= R_max[bus,1]);
#         # Initial battery SOC
#         @constraint(m, B_rt[bus,1]==B_feedback[bus,1]);
#
#         # SOC constrains on real and reactive power on bus
#         maxS_rt = sqrt(1+Q_gamma^2)*(Pd[bus,1]-R_max[bus,1]-
#         positive_scalar(
#         icdf*sqrt(pd.sigma[bus,1]+pg.sigma[bus,1])+pg.mu[bus,1]));
#         if maxS_rt>=0
#             @constraint(m,
#             [maxS_rt, 0.5*maxS_rt,
#             P_bus_rt[bus,1], Q_bus_rt[bus,1]] in RotatedSecondOrderCone())
#         else
#             @constraint(m,
#             [-maxS_rt, -0.5*maxS_rt,
#             P_bus_rt[bus,1], Q_bus_rt[bus,1]] in RotatedSecondOrderCone())
#         end
#         # a list of branchs that FROM the bus of interest
#         sub_branch_list = findall(one->one==1, C_ind[bus,:]) # Out
#         # a list of branchs that POINT TO the bus of interest
#         add_branch_list = findall(minusone->minusone==-1, C_ind[bus,:]) #In
#
#         # real power on bus is demand minus solar and discharge
#         # AKA power injection
#         @constraint(m, P_bus_rt[bus,1]==
#         Pd[bus,1]
#         -(Pg_rt[bus,1]+R_rt[bus,1]));
#
#         if bus_struct.type[bus]==1 # non-generator bus
#             # box constraint on voltage
#             @constraint(m, V_min^2<= v_rt[bus,1]);
#             @constraint(m, v_rt[bus,1]<= V_max^2);
#
#
#
#             # Power Balance Equations
#             # Power injection = Power Flow In - Power Flow Out
#             # if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
#             #     @constraint(m, P_bus_rt[bus,1]==
#             #         sum(P_br_rt[branch,1] for branch in add_branch_list)
#             #         -sum(P_br_rt[branch,1] for branch in sub_branch_list)
#             #         );
#             #     @constraint(m, Q_bus_rt[bus,1]==
#             #         sum(Q_br_rt[branch,1] for branch in add_branch_list)
#             #         -sum(Q_br_rt[branch,1] for branch in sub_branch_list));
#             # elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
#             #
#             #     @constraint(m, P_bus_rt[bus,1]==
#             #         sum(P_br_rt[branch,1] for branch in add_branch_list)
#             #         );
#             #     @constraint(m, Q_bus_rt[bus,1]==
#             #         sum(Q_br_rt[branch,1] for branch in add_branch_list));
#             # elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
#             #     @constraint(m, P_bus_rt[bus,1]==
#             #         -sum(P_br_rt[branch,1] for branch in sub_branch_list)
#             #         );
#             #     @constraint(m, Q_bus_rt[bus,1]==
#             #         -sum(Q_br_rt[branch,1] for branch in sub_branch_list));
#             # end
#         else # generator bus
#             ########
#             @constraint(m, v_rt[bus,1]== 1);
#
#
#             # only one generator, and is regraded as the slack bus
#             ########
#             # @constraint(m, V_min^2<= v_rt[bus,1]);
#             # @constraint(m, v_rt[bus,1]<= V_max^2);
#
#             # identify the id of the generator that connects to the bus of interest
#             gen_id = findall(bus_id ->bus_id == bus, gen_struct.bus);
#             # Power Balance Equations
#             # Power injection = Power Flow In - Power Flow Out + Power Generatered
#             # if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
#             #     @constraint(m, P_bus_rt[bus,1]==
#             #         sum(P_br_rt[branch,1] for branch in add_branch_list)
#             #         -sum(P_br_rt[branch,1] for branch in sub_branch_list)
#             #         +P_gen_rt[gen_id[1], 1]);
#             #     @constraint(m, Q_bus_rt[bus,1]==
#             #         sum(Q_br_rt[branch,1] for branch in add_branch_list)
#             #         -sum(Q_br_rt[branch,1] for branch in sub_branch_list)
#             #         +Q_gen_rt[gen_id[1], 1]);
#             # elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
#             #     @constraint(m, P_bus_rt[bus,1]==
#             #         sum(P_br_rt[branch,1] for branch in add_branch_list)
#             #         +P_gen_rt[gen_id[1], 1]);
#             #     @constraint(m, Q_bus_rt[bus,1]==
#             #         sum(Q_br_rt[branch,1] for branch in add_branch_list)
#             #         +Q_gen_rt[gen_id[1], 1]);
#             # elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
#             #     @constraint(m, P_bus_rt[bus,1]==
#             #         -sum(P_br_rt[branch,1] for branch in sub_branch_list)
#             #         +P_gen_rt[gen_id[1], 1]);
#             #     @constraint(m, Q_bus_rt[bus,1]==
#             #         -sum(Q_br_rt[branch,1] for branch in sub_branch_list)
#             #         +Q_gen_rt[gen_id[1], 1]);
#             # end
#         end
#     end
#     # voltage constraint for LinDistFlow
#     for branch =1:NoBr
#         # println(string("branch # ", branch))
#         fbus = branch_struct.fbus[branch];
#         tbus = branch_struct.tbus[branch];
#         tbus_branch_list = findall(isafbus->isafbus==tbus, branch_struct.fbus[:,1]);
#         # println("tbus branch list:")
#         # println(tbus_branch_list)
#         # @constraint(m, v_rt[fbus]-v_rt[tbus]==
#         #     2*(r[branch]*P_br_rt[branch]+x[branch]*Q_br_rt[branch]))
#         if isempty(tbus_branch_list)
#             @constraint(m,P_br_rt[branch,1]==
#             r[branch]*l_br_rt[branch,1]+P_bus_rt[tbus, 1]);
#
#             @constraint(m,Q_br_rt[branch,1]==
#             x[branch]*l_br_rt[branch,1]+Q_bus_rt[tbus, 1]);
#         else
#             @constraint(m, P_br_rt[branch,1]==
#             r[branch]*l_br_rt[branch,1]+P_bus_rt[tbus, 1]
#             +sum(P_br_rt[branch,1] for branch in tbus_branch_list));
#
#             @constraint(m, Q_br_rt[branch,1]==
#             x[branch]*l_br_rt[branch,1]+Q_bus_rt[tbus, 1]
#             +sum(Q_br_rt[branch,1] for branch in tbus_branch_list));
#         end
#
#         @constraint(m,v_rt[tbus,1] ==
#         v_rt[fbus,1] + (r[branch]^2+x[branch]^2)*l_br_rt[branch,1]
#         -2*(r[branch]*P_br_rt[branch,1]+x[branch]*Q_br_rt[branch,1]));
#
#         @constraint(m, [0.5*l_br_rt[branch,1], v_rt[fbus,1],
#         P_br_rt[branch,1], Q_br_rt[branch,1]] in RotatedSecondOrderCone());
#     end
#
#     # box constraint on generator
#     for Gen=1:NoGen
#         tbus_branch_list = findall(isafbus->isafbus==1, branch_struct.fbus[:,1]);
#         # println("generator tbus branch list:")
#         # println(tbus_branch_list)
#         # println(tbus_branch_list)
#         @constraint(m,
#             P_gen_rt[Gen,1] == P_bus_rt[1,1]
#             +sum(P_br_rt[branch,1] for branch in tbus_branch_list));
#         @constraint(m,
#             Q_gen_rt[Gen,1] == Q_bus_rt[1,1]
#             +sum(Q_br_rt[branch,1] for branch in tbus_branch_list));
#         @constraint(m,
#             P_gen_rt[Gen,1]<=gen_struct.Pmax[Gen]/baseMVA)
#         @constraint(m,
#             P_gen_rt[Gen,1]>=gen_struct.Pmin[Gen]/baseMVA)
#         @constraint(m,
#             Q_gen_rt[Gen,1]<=gen_struct.Qmax[Gen]/baseMVA)
#         @constraint(m,
#             Q_gen_rt[Gen,1]>=gen_struct.Qmin[Gen]/baseMVA)
#     end
#
#     # bus
#     @variable(m, Pg[1:SN, 1:NoBus, 1:T-1]); # the real power output
#     @variable(m, B[1:SN, 1:NoBus, 1:T-1]); # the storage
#     @variable(m, R[1:SN, 1:NoBus, 1:T-1]);# the charge/discharge rate
#     @variable(m, P_bus[1:SN, 1:NoBus, 1:T-1])
#     @variable(m, Q_bus[1:SN, 1:NoBus, 1:T-1])
#     @variable(m, v[1:SN, 1:NoBus, 1:T-1])
#
#     # Branch level
#     @variable(m, P_br[1:SN,1:NoBr, 1:T-1])
#     @variable(m, Q_br[1:SN,1:NoBr, 1:T-1])
#     @variable(m, l_br[1:SN,1:NoBr, 1:T-1])
#     # Gen Level
#     @variable(m, P_gen[1:SN, 1:NoGen, 1:T-1])
#     @variable(m, Q_gen[1:SN, 1:NoGen, 1:T-1])
#
#     # for different scenario
#     for scenario = 1:SN
#         # for different time at prediction horizion
#         for t=1:T-1
#             # for different bus
#             for bus=1:NoBus
#                 # box constraints on Solar (Pg), Battery discharge (R)
#                 @constraint(m, Pg_min[bus, t]<=Pg[scenario, bus ,t]);
#                 @constraint(m, Pg[scenario, bus,t]<=
#                     positive_scalar(
#                     icdf*sqrt(pd.sigma[bus,t+1]+pg.sigma[bus,t+1])+pg.mu[bus,t+1])
#                     );
#                 @constraint(m, R_min[bus,t+1] <= R[scenario,bus,t]);
#                 @constraint(m, R[scenario,bus,t]<= R_max[bus,t+1]);
#                 # battery box constraint on SOC
#                 @constraint(m, B_min[bus,t+1] <= B[scenario,bus,t]);
#                 @constraint(m, B[scenario,bus,t]<= B_max[bus,t+1]);
#                 # battery SOC dynamics
#                 if t==1
#                     @constraint(m, B[scenario,bus,1] == B_rt[bus,1]
#                         -delta_t*(R_rt[bus,1]*baseMVA))
#                 else
#                     @constraint(m, B[scenario,bus,t] ==
#                         B[scenario,bus,t-1] - R[scenario,bus,t-1]*baseMVA*delta_t)
#                 end
#                 # real power on bus is demand minus solar and discharge
#                 # AKA power injection
#                 @constraint(m, P_bus[scenario,bus,t]==
#                 Pd[bus,t+1]
#                 -Pg[scenario,bus,t]-R[scenario, bus,t]);
#
#                 # SOC constrains on real and reactive power on bus
#                 # maxS = sqrt(1+Q_gamma^2)*(Pd[bus,t+1]-positive_scalar(
#                 # icdf*sqrt(pd.sigma[bus,t+1]+pg.sigma[bus,t+1])+pg.mu[bus,t+1])
#                 # -R_max[bus,t+1]);
#                 maxS = sqrt(1+Q_gamma^2)*(Pd[bus,t+1]-positive_scalar(
#                 icdf*sqrt(pd.sigma[bus,t+1]+pg.sigma[bus,t+1])+pg.mu[bus,t+1])
#                 -R_max[bus,t+1]);
#                 if maxS >=0
#                     @constraint(m,
#                     [maxS, 0.5*maxS, P_bus[scenario,bus,t], Q_bus[scenario,bus,t]]
#                      in RotatedSecondOrderCone())
#                 else
#                     @constraint(m,
#                     [-maxS, -0.5*maxS, P_bus[scenario,bus,t], Q_bus[scenario,bus,t]]
#                      in RotatedSecondOrderCone())
#                 end
#
#                 # a list of branchs that FROM the bus of interest
#                 sub_branch_list = findall(one->one==1, C_ind[bus,:]) # out
#                 # a list of branchs that POINT TO the bus of interest
#                 add_branch_list = findall(minusone->minusone==-1, C_ind[bus,:]) # In
#
#                 if bus_struct.type[bus]==1 # non-generator bus
#                     # box constraint on bus
#                     @constraint(m, V_min^2<= v[scenario,bus,t]);
#                     @constraint(m, v[scenario,bus,t]<= V_max^2);
#                     # Power Balance Equations
#                     # Power injection = Power Flow In - Power Flow Out
#                     # if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
#                     #     @constraint(m, P_bus[scenario,bus,t]==
#                     #         sum(P_br[scenario,branch,t] for branch in add_branch_list)
#                     #         -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
#                     #         );
#                     #     @constraint(m, Q_bus[scenario,bus,t]==
#                     #         sum(Q_br[scenario,branch,t] for branch in add_branch_list)
#                     #         -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
#                     #         );
#                     # elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
#                     #     @constraint(m, P_bus[scenario,bus,t]==
#                     #         sum(P_br[scenario,branch,t] for branch in add_branch_list)
#                     #         );
#                     #     @constraint(m, Q_bus[scenario,bus,t]==
#                     #         sum(Q_br[scenario,branch,t] for branch in add_branch_list)
#                     #         );
#                     # elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
#                     #
#                     #     @constraint(m, P_bus[scenario,bus,t]==
#                     #         -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
#                     #         );
#                     #     @constraint(m, Q_bus[scenario,bus,t]==
#                     #         -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
#                     #         );
#                     # end
#                 else # genertor bus
#
#                     # @constraint(m, V_min^2<= v[scenario,bus,t]);
#                     # @constraint(m, v[scenario,bus,t]<= V_max^2);
#                     gen_id = findall(bus_id ->bus_id == bus, gen_struct.bus);
#                     @constraint(m, v[scenario,bus,t]==1);
#                     # if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
#                     #     @constraint(m, P_bus[scenario,bus,t]==
#                     #         sum(P_br[scenario,branch,t] for branch in add_branch_list)
#                     #         -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
#                     #         +P_gen[scenario,gen_id[1],t]);
#                     #     @constraint(m, Q_bus[scenario,bus,t]==
#                     #         sum(Q_br[scenario,branch,t] for branch in add_branch_list)
#                     #         -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
#                     #         +Q_gen[scenario,gen_id[1],t]);
#                     # elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
#                     #     @constraint(m, P_bus[scenario,bus,t]==
#                     #         sum(P_br[scenario,branch,t] for branch in add_branch_list)
#                     #         +P_gen[scenario,gen_id[1],t]);
#                     #     @constraint(m, Q_bus[scenario,bus,t]==
#                     #         sum(Q_br[scenario,branch,t] for branch in add_branch_list)
#                     #         +Q_gen[scenario,gen_id[1],t]);
#                     # elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
#                     #     @constraint(m, P_bus[scenario,bus,t]==
#                     #         -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
#                     #         +P_gen[scenario,gen_id[1],t]);
#                     #     @constraint(m, Q_bus[scenario,bus,t]==
#                     #         -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
#                     #         +Q_gen[scenario,gen_id[1],t]);
#                     # end
#                 end
#             end
#
#             for branch =1:NoBr
#                 fbus = Int(branch_struct.fbus[branch]);
#                 tbus = Int(branch_struct.tbus[branch]);
#                 tbus_branch_list = findall(isafbus->isafbus==tbus, branch_struct.fbus[:,1]);
#                 # @constraint(m, v[scenario,fbus,t]-v[scenario,tbus,t]==
#                 #     2*(r[branch]*P_br[scenario,branch,t]
#                 #     +x[branch]*Q_br[scenario,branch,t]))
#                 if isempty(tbus_branch_list)
#                     # println(branch)
#                     @constraint(m,P_br[scenario,branch,t]==
#                     r[branch]*l_br[scenario,branch,t]+P_bus[scenario, tbus, t]);
#
#                     @constraint(m,Q_br[scenario,branch,t]==
#                     x[branch]*l_br[scenario,branch,t]+Q_bus[scenario,tbus, t]);
#                 else
#                     @constraint(m,P_br[scenario,branch,t]==
#                     r[branch]*l_br[scenario,branch,t]+P_bus[scenario, tbus, t]
#                     +sum(P_br[scenario,branch,t] for branch in tbus_branch_list));
#
#                     @constraint(m,Q_br[scenario,branch,t]==
#                     x[branch]*l_br[scenario,branch,t]+Q_bus[scenario,tbus, t]
#                     +sum(Q_br[scenario,branch,t] for branch in tbus_branch_list));
#                 end
#
#                 @constraint(m,v[scenario,tbus,t]==
#                 v[scenario,fbus,t]+(r[branch]^2+x[branch]^2)*l_br[scenario,branch,t]-
#                 2*(r[branch]*P_br[scenario,branch,t]+x[branch]*Q_br[scenario,branch,t]));
#
#                 @constraint(m, [0.5*l_br[scenario,branch,t], v[scenario,fbus,t],
#                 P_br[scenario,branch,t], Q_br[scenario,branch,t]] in RotatedSecondOrderCone());
#             end
#
#             for Gen=1:NoGen
#                 tbus_branch_list = findall(isafbus->isafbus==1, branch_struct.fbus[:,1]);
#                 @constraint(m,
#                     P_gen[scenario,Gen,t] == P_bus[scenario,1,t]+
#                     sum(P_br[scenario,branch,t] for branch in tbus_branch_list));
#                 @constraint(m,
#                     Q_gen[scenario,Gen,t] == Q_bus[scenario,1,t]+
#                     sum(Q_br[scenario,branch,t] for branch in tbus_branch_list));
#                 @constraint(m,
#                     P_gen[scenario,Gen,t]<=gen_struct.Pmax[Gen]/baseMVA)
#                 @constraint(m,
#                     P_gen[scenario,Gen,t]>=gen_struct.Pmin[Gen]/baseMVA)
#                 @constraint(m,
#                     Q_gen[scenario,Gen,t]<=gen_struct.Qmax[Gen]/baseMVA)
#                 @constraint(m,
#                     Q_gen[scenario,Gen,t]>=gen_struct.Qmin[Gen]/baseMVA)
#                 # @constraint(m,
#                 #     Q_gen[scenario,Gen,t]<=1)
#                 # @constraint(m,
#                 #     Q_gen[scenario,Gen,t]>=-1)
#             end
#         end
#     end
#
#
#     if ancillary_type == "10min" || ancillary_type == "30min"
#             # println(Real_Bus_list)
#             # RT
#             @variable(m, P_rsrv_rt)
#             @variable(m, B_rsrv_rt)
#
#             @constraint(m, P_rsrv_rt>=0)
#             @constraint(m, B_rsrv_rt>=0)
#
#             if current_time <= tau
#                 @constraint(m, B_rsrv_rt==0)
#
#             elseif current_time - tau>=1
#                 ini_fb = max(current_time-tau-k+1,1);
#                 fni_fb = current_time-tau;
#                 length_fb = fni_fb - ini_fb +1;
#                 mult_fb = k-length_fb+1:k;
#                 # println(ini_fb)
#                 # println
#                 temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb]
#                 if length_fb >= 2
#                     for f_rsrv_fb_n=2:length_fb;
#                         temp_f_rsrv_c_fb = (temp_f_rsrv_c_fb+
#                             mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n]);
#                     end
#                 end
#                 @constraint(m, B_rsrv_rt==floor(delta_t*temp_f_rsrv_c_fb*1000)/1000)
#                 # println("nominal B_rsrv")
#                 # println(floor(delta_t*temp_f_rsrv_c_fb*1000)/1000)
#                 # println("sum B_feedback")
#                 # println(sum(B_feedback))
#             end
#             # for Bus in Real_Bus_list
#             #     list = setdiff(Real_Bus_list, Bus)
#             #     # println(list)
#             #     # println(sum(B_feedback[bus, 1] for bus in list))
#             #     @constraint(m, B_rsrv_rt <= sum(B_rt[bus, 1] for bus in list))
#             # end
#
#             @constraint(m, B_rsrv_rt <= sum(B_rt))
#             # println(sum(B_feedback[bus, 1] for bus in Non_affected_Bus_list))
#             # @constraint(m, B_rsrv_rt <=
#             #     sum(B_rt[bus, 1] for bus in Non_affected_Bus_list))
#             # @constraint(m, B_rsrv_rt <= 4)
#             # Scenario
#             @variable(m, P_rsrv[1:SN,1:T-1])
#             @variable(m, B_rsrv[1:SN,1:T-1])
#             # @constraint(m, P_rsrv==3)
#
#             for scenario=1:SN
#                 for t_real=current_time+1:current_time+T-1
#                     @constraint(m, P_rsrv[scenario, t_real-current_time]>=0)
#                     @constraint(m, B_rsrv[scenario, t_real-current_time]>=0)
#                     ini = t_real-tau-k+1;
#                     fin = t_real-tau;
#                     if fin <=0
#                         @constraint(m, B_rsrv[scenario, t_real-current_time]==0)
#                     elseif fin < current_time && fin>=1
#                         ini_fb = max(1,ini);
#                         fin_fb = fin;
#                         length_fb = fin_fb-ini_fb+1;
#                         mult_fb = k-length_fb+1:k;
#                         temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb];
#                         if length_fb >= 2
#                             for f_rsrv_fb_n=2:length_fb;
#                                 temp_f_rsrv_c_fb = temp_f_rsrv_c_fb+
#                                     mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n];
#                             end
#                         end
#                         @constraint(m, B_rsrv[scenario, t_real-current_time]==
#                             floor(delta_t*temp_f_rsrv_c_fb*1000)/1000);
#                     elseif fin == current_time
#                         if current_time == 1
#                             @constraint(m, B_rsrv[scenario, t_real-current_time]==
#                             delta_t*k*P_rsrv_rt);
#                         else
#                             ini_fb = max(1, ini);
#                             fin_fb = fin-1;
#                             length_fb = fin_fb-ini_fb+1;
#                             mult_fb = k-length_fb:k-1;
#                             temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb];
#                             if length_fb >= 2
#                                 for f_rsrv_fb_n=2:length_fb;
#                                     temp_f_rsrv_c_fb = temp_f_rsrv_c_fb+
#                                        mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n];
#                                  end
#                             end
#                             # @constraint(m, B_rsrv[scenario, t_real-current_time]==
#                             #     floor(delta_t*k*P_rsrv_rt*1000)/1000);
#                             @constraint(m, B_rsrv[scenario, t_real-current_time]==
#                                  delta_t*k*P_rsrv_rt
#                                  +floor(delta_t*temp_f_rsrv_c_fb*1000)/1000);
#                         end
#                     elseif ini < current_time && fin > current_time
#                         if current_time == 1
#                             ini_sc = current_time+1;
#                             fin_sc = fin;
#                             length_sc = fin_sc - ini_sc +1;
#                             mult_sc = k-length_sc+1:k;
#                             mult_rt = k-length_sc;
#                              temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time]
#                             if length_sc > 1
#                                 for f_rsrv_sc_n=2:length_sc
#                                     temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+
#                                         mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
#                                 end
#                             end
#                             @constraint(m, B_rsrv[scenario, t_real-current_time] ==
#                                 delta_t*mult_rt*P_rsrv_rt
#                                 +delta_t*temp_f_rsrv_c_sc);
#                         else
#                             ini_sc = current_time+1;
#                             fin_sc = fin;
#                             length_sc = fin_sc - ini_sc +1;
#                             mult_sc = k-length_sc+1:k;
#                             temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time]
#                             if length_sc > 1
#                                 for f_rsrv_sc_n=2:length_sc
#                                     temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+
#                                         mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
#                                 end
#                             end
#                             mult_rt = k-length_sc;
#                             ini_fb = max(1, ini);
#                             fin_fb = current_time-1;
#                             length_fb = fin_fb-ini_fb+1;
#                             temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb];
#                             if length_fb >= 2
#                                 for f_rsrv_fb_n=2:length_fb;
#                                     temp_f_rsrv_c_fb = temp_f_rsrv_c_fb+
#                                         mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n];
#                                 end
#                             end
#                             @constraint(m, B_rsrv[scenario, t_real-current_time] ==
#                             delta_t*mult_rt*P_rsrv_rt+
#                             delta_t*temp_f_rsrv_c_sc+
#                             floor(delta_t*temp_f_rsrv_c_fb*1000)/1000);
#                         end
#                     elseif ini == current_time
#                         ini_sc = current_time+1;
#                         fin_sc = fin;
#                         mult_sc = 2:k;
#                         temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time];
#                         for f_rsrv_sc_n=2:k-1
#                             temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+
#                                 mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
#                         end
#                         @constraint(m, B_rsrv[scenario, t_real-current_time] ==
#                             delta_t*P_rsrv_rt
#                             +delta_t*temp_f_rsrv_c_sc);
#                     elseif ini > current_time
#                         ini_sc = ini;
#                         fin_sc = fin;
#                         mult_sc = 1:k;
#                         temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time]
#                         for f_rsrv_sc_n=2:k
#                             temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+
#                                 mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
#                         end
#                         @constraint(m, B_rsrv[scenario, t_real-current_time]
#                             == delta_t*temp_f_rsrv_c_sc);
#
#                         if  t_real-current_time >= T-tau-1
#                             @constraint(m, P_rsrv[scenario, t_real-current_time]==0)
#                         end
#                     end
#                     # for Bus in Real_Bus_list
#                     #     list = setdiff(Real_Bus_list, Bus)
#                     #     @constraint(m, B_rsrv[scenario, t_real-current_time]
#                     #         <= sum(B[scenario, bus, t_real-current_time] for bus in list))
#                     # end
#
#                     @constraint(m, B_rsrv[scenario, t_real-current_time]
#                          <= sum(B[scenario, :, t_real-current_time]))
#
#                 end
#             end
#     end
#     # #
#
#     if ancillary_type == "10min" || ancillary_type == "30min"
#         @objective(m, Min,
#             fn_cost_RHC_anc(delta_t, P_gen, P_gen_rt, Pg_rt,Pg, P_rsrv_rt,
#             P_rsrv,price,pg,pd,beta,SN,obj, gen_struct, bus_struct, branch_struct))
#     else
#         @objective(m, Min,
#             fn_cost_RHC_rt(delta_t, P_gen, P_gen_rt, Pg_rt,Pg, price,
#                 pg, pd, beta, SN, obj, gen_struct, bus_struct, branch_struct, baseMVA))
#     end
#     #
#     status=optimize!(m);
#     println(string("    ----", termination_status(m)))
#     terminate_s = termination_status(m);
#     # println(MOI.PrimalStatus())
#     # println(MOI.DualStatus())
#     cost_o = JuMP.objective_value(m);
#     time_solve=MOI.get(m, MOI.SolveTime());
#     println(string("    ----Solve Time: ", time_solve))
#     println(string("    ----Optimal Cost of the RHC: ", cost_o))
#
#     println(string("    ----aggregate demand: ", baseMVA*sum(Pd[:,1])))
#     ###############
#     R_rt_o=JuMP.value.(R_rt)
#     R_o=JuMP.value.(R)
#     R_sum_ct = sum(R_rt_o[:, 1]*baseMVA)
#     println(string("    ----aggregate battery disharge: ", R_sum_ct))
#     R_traj = hcat(R_rt_o, R_o[1,:,:])
#     ###############
#     P_br_rt_o=JuMP.value.(P_br_rt)
#     P_br_o=JuMP.value.(P_br)
#     P_br_traj = hcat(P_br_rt_o, P_br_o[1,:,:])
#
#     ###############
#     Q_br_rt_o=JuMP.value.(Q_br_rt)
#     Q_br_o=JuMP.value.(Q_br)
#     Q_br_traj = hcat(Q_br_rt_o, Q_br_o[1,:,:])
#     ###############
#     v_rt_o=JuMP.value.(v_rt)
#     v_o=JuMP.value.(v)
#     v_traj = hcat(v_rt_o, v_o[1,:,:])
#
#     ################
#     B_rt_o=JuMP.value.(B_rt)
#     B_o=JuMP.value.(B)
#     B_traj = hcat(B_rt_o, B_o[1,:,:])
#     ###############
#     Pg_rt_o=JuMP.value.(Pg_rt)
#     Pg_o=JuMP.value.(Pg)
#     Pg_traj = hcat(Pg_rt_o, Pg_o[1,:,:])
#
#     ############
#     P_bus_rt_o=JuMP.value.(P_bus_rt)
#     P_bus_o=JuMP.value.(P_bus)
#     P_bus_traj = hcat(P_bus_rt_o, P_bus_o[1,:,:])
#
#     ############
#     Q_bus_rt_o=JuMP.value.(Q_bus_rt)
#     Q_bus_o=JuMP.value.(Q_bus)
#     Q_bus_traj = hcat(Q_bus_rt_o, Q_bus_o[1,:,:])
#     println(string("    ----aggregate reactive on bus: ",
#         sum(Q_bus_traj)))
#
#     ############
#     P_gen_rt_o=JuMP.value.(P_gen_rt)
#     P_gen_o=JuMP.value.(P_gen)
#     P_gen_traj = hcat(P_gen_rt_o, P_gen_o[1,:,:])
#     P_gen_sum_traj = hcat(sum(P_gen_rt_o), sum(P_gen_o[1,:,:], dims=1))
#     P_gen_sum_ct = sum(P_gen_rt_o[gen, 1]*baseMVA for gen in 1:NoGen)
#     println(string("    ----current real generation: ", P_gen_sum_ct))
#     println(string("    ----aggregate real generation: ", sum(P_gen_sum_traj)))
#
#     ###############
#     l_br_rt_o=JuMP.value.(l_br_rt)
#     l_br_o=JuMP.value.(l_br)
#     l_br_traj = hcat(l_br_rt_o, l_br_o[1,:,:])
#
#     ############
#     Q_gen_rt_o = JuMP.value.(Q_gen_rt)
#     Q_gen_o = JuMP.value.(Q_gen)
#     Q_gen_traj = hcat(Q_gen_rt_o, Q_gen_o[1,:,:])
#     Q_gen_sum_traj = hcat(sum(Q_gen_rt_o), sum(Q_gen_o[1,:,:], dims=1))
#     Q_gen_sum_ct = sum(Q_gen_rt_o[gen, 1]*baseMVA for gen in 1:NoGen)
#     println(string("    ----current reactive generation: ", Q_gen_sum_ct))
#     println(string("    ----aggregate reactive generation: ",
#         sum(Q_gen_sum_traj)))
#     ###########
#
#
#     # println(sum(P_bus_rt_o))
#     # println(sum(Pd[:,1]))
#     # println(sum(P_bus_o[1,:,1]))
#     # println(sum(Pd[:,2]))
#     # println(P_gen_rt_o[1,1])
#     # println(sum(Pd[:,1]))
#     # println(P_gen_o[1,1,1])
#     # println(sum(Pd[:,2]))
#     # for branch =1:NoBr
#     #     println(string("branch",branch))
#     #     t=2;
#     #     scenario=1;
#     #     fbus = Int(branch_struct.fbus[branch]);
#     #     tbus = Int(branch_struct.tbus[branch]);
#     #     tbus_branch_list = findall(isafbus->isafbus==tbus, branch_struct.fbus[:,1]);
#     #     if isempty(tbus_branch_list)
#     #         println(P_br_o[scenario,branch,t]-
#     #         r[branch]*l_br_o[scenario,branch,t]-P_bus_o[scenario, tbus, t])
#     #         println(P_br_o[scenario,branch,t])
#     #         println(r[branch]*l_br_o[scenario,branch,t])
#     #         println(P_bus_o[scenario, tbus, t])
#     #     else
#     #         println(P_br_o[scenario,branch,t]-
#     #         r[branch]*l_br_o[scenario,branch,t]-P_bus_o[scenario, tbus, t]
#     #         -sum(P_br_o[scenario,branch,t] for branch in tbus_branch_list))
#     #         println(P_br_o[scenario,branch,t])
#     #         println(r[branch]*l_br_o[scenario,branch,t])
#     #         println(P_bus_o[scenario, tbus, t])
#     #         println(sum(P_br_o[scenario,branch,t] for branch in tbus_branch_list))
#     #     end
#     # end
#     # println(sum(P_bus_o[1,:,1]))
#     # println(sum(Pd[:,2]))
#     # println("print voltage constraint")
#     # for branch = 17:17
#     #     fbus = branch_struct.fbus[branch];
#     #     tbus = branch_struct.tbus[branch];
#     #     println(Q_br_rt_o[branch,1]-
#     #     x[branch]*l_br_rt_o[branch,1]-Q_bus_rt_o[tbus, 1])
#     #     # println(v_rt_o[fbus]-v_rt_o[tbus]-
#     #     #     2*(r[branch]*P_br_rt_o[branch,1]))
#     # end
#     # for t=1:287
#     #     for branch = 17:17
#     #         fbus = Int(branch_struct.fbus[branch]);
#     #         tbus = Int(branch_struct.tbus[branch]);
#     #         # println(Q_bus_o[1,tbus, t])
#     #         println(Q_br_o[1,branch, t]-
#     #         x[branch]*l_br_o[1,branch, t]-Q_bus_o[1,tbus, t])
#     #         # println(v_o[1,fbus,t]-v_o[1,tbus,t]-
#     #         #     2*(r[branch]*P_br_o[1,branch,t]
#     #         #     ))
#     #     end
#     # end
#
#     ############
#     if ancillary_type == "10min" || ancillary_type == "30min"
#         P_rsrv_rt_o=JuMP.value(P_rsrv_rt);
#         P_rsrv_s=JuMP.value.(P_rsrv);
#         P_rsrv_total = hcat(P_rsrv_rt_o[1,1], reshape(P_rsrv_s[1,:], 1, 287));
#         B_rsrv_rt_o=JuMP.value(B_rsrv_rt);
#         B_rsrv_s=JuMP.value.(B_rsrv);
#         B_rsrv_total = hcat(B_rsrv_rt_o[1,1], reshape(B_rsrv_s[1,:], 1, 287));
#     else
#         P_rsrv_rt_o = 0;
#         P_rsrv_total = zeros(1,T);
#         B_rsrv_rt_o = 0;
#         B_rsrv_total = zeros(1,T);
#     end
#     # println(P_rsrv_rt_o)
#     # println
#     if sum(Pg_rt_o)<=sum(pg.mu_ct)
#     ############
#         Cost_real = delta_t*(sum(P_gen_rt_o)*price.lambda_ct
#             +beta*(sum(Pg_rt_o)-sum(pg.mu_ct)))*baseMVA-
#             delta_t*P_rsrv_rt_o*price.alpha_ct;
#     else
#         Cost_real = delta_t*(sum(P_gen_rt_o)*price.lambda_ct
#             +price.lambda_ct*(sum(Pg_rt_o)-sum(pg.mu_ct)))*baseMVA-
#             delta_t*P_rsrv_rt_o*price.alpha_ct;
#     end
#    # println(price.alpha_ct)
#     P_cul = sum(pg.mu_ct)-sum(Pg_rt_o)
#     println(string("    ----curtailment solar: ", P_cul))
#     alpha_1 = price.alpha_scenario[1,:]
#     lambda_1 = price.lambda_scenario[1,:]
#
#     println(string("    ----Optimzal cost at this instance: ", Cost_real))
#     pg_upper = hcat(
#     icdf*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])+pg.mu[:,:])
#
#     val_opt = (R=(R_rt_o), B=(B_rt_o), Pg=(Pg_rt_o), P_gen=(P_gen_rt_o), Q_gen=(Q_gen_rt_o),
#         P_rsrv=(P_rsrv_rt_o), B_rsrv=(B_rsrv_rt_o), P_bus=(P_bus_rt_o), Q_bus=(Q_bus_rt_o),
#         v=(v_rt_o), P_br=(P_br_rt_o), Q_br=(Q_br_rt_o),
#         Cost_real=(Cost_real), time_solve=(time_solve),
#         P_rsrv_total=(P_rsrv_total), B_rsrv_total=(B_rsrv_total), alpha_1=(alpha_1),
#         R_traj = (R_traj), B_traj=(B_traj), Pg_traj=(Pg_traj), P_gen_traj=(P_gen_traj),
#         Q_gen_traj=(Q_gen_traj),
#         P0_traj = (P_gen_traj), P_bus_traj=(P_bus_traj), Q_bus_traj=(Q_bus_traj),
#         pg_upper=(pg_upper),
#         lambda_1=(lambda_1), pd=(pd), pg=(pg), v_traj=(v_traj), P_br_traj=(P_br_traj),
#         Q_br_traj=(Q_br_traj), terminate_s = (terminate_s), P_cul = (P_cul))
#
#
#     return val_opt
# end
