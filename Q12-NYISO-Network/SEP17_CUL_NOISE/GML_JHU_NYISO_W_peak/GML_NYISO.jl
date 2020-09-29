function GML_solver(date, current_time, B_inp, raw_data)

    # sys para
    # ancillary_type="10min"
    # ancillary_type="30min"
    ancillary_type="without";
    # # of timeslots
    T=288;
    # # of secnarios
    SN=3;
    # solar peneration rate
    p_rate = 0.25;
    ct_printout = string("===== Solar penetration rate", p_rate);
    println("=================================================")
    println(ct_printout)
    ################################################################################
    # penetation level
    icdf = 0;
    Pred_length = 12;


    solar_error_max = 0.025;
    # #
    if solar_error_max == 0.025
        pg_noise = matread("GML_data/price-solar-demand/solar_noise_0025.mat")["solar_noise"];
    elseif solar_error_max == 0.05
        pg_noise = matread("GML_data/price-solar-demand/solar_noise_005.mat")["solar_noise"];
    elseif solar_error_max == 0.1
        pg_noise = matread("GML_data/price-solar-demand/solar_noise_01.mat")["solar_noise"];
    elseif solar_error_max == 0
        pg_noise = zeros(1,12);
    end

    shunt_struct = raw_data.shunt_struct;
    bus_struct = raw_data.bus_struct;
    branch_struct = raw_data.branch_struct;
    gen_struct = raw_data.gen_struct;
    price_raw = raw_data.price_raw;
    delta_rt_raw = raw_data.delta_rt_raw;
    pd_raw = raw_data.pd_raw;
    pd_noise  = raw_data.pd_noise;
    pg_raw = raw_data.pg_raw;

    baseMVA = 100;

    NoShunt = length(shunt_struct.P)
    NoBus = length(bus_struct.baseKV)


    B_cap = 75;

    if current_time == 1
        P_rsrv_feedback = [];
        B_feedback=zeros(NoShunt, 1);
        B_feedback[2,1] = B_inp[2,1];
        B_feedback[4,1] = B_inp[1,1];
        B_feedback[5,1] = B_inp[3,1];
    else
        B_feedback = read_B_out(current_time-1)
        # B_feedback=zeros(NoShunt, 1);
        # B_feedback[2,1] = B_inp[2,1];
        # B_feedback[4,1] = B_inp[1,1];
        # B_feedback[5,1] = B_inp[3,1];
        if ancillary_type == "10min" || ancillary_type == "30min"
            P_rsrv_feedback = read_RSRV_out()
            # P_rsrv_feedback = zeros(1, current_time-1)
        else
            P_rsrv_feedback = zeros(1, current_time-1)
        end
    end

    feedback = (B_feedback=(B_feedback), P_rsrv_feedback=(P_rsrv_feedback));




    # price = price_traj_det(current_time, ancillary_type, price_raw, T);
    price = price_traj(current_time, ancillary_type, price_raw, delta_rt_raw, T, Pred_length)

    # println(price.lambda_scenario[1,1])
    # ################################################################################
    # # plotting prices
    # ################################################################################

    # plot(1:288, reshape(price.lambda_scenario[1,:], 288,1), label="lambda")
    # plot(1:288, reshape(price.alpha_scenario[1,:], 288,1), label="alpha")
    # # # lambda = reshape(price.lambda_scenario[1,:], 288,1)
    # # # alpha = reshape(price.alpha_scenario[1,:], 288,1)
    # # # price = hcat(lambda, alpha)
    # # # CSV.write("price.csv", DataFrame(price, [:lambda, :alpha]))
    # # # println(sum(price.alpha_scenario[1,:]))
    # ################################################################################
    # pd = pd_traj_pu_det(current_time, pd_raw, T,  NoShunt, baseMVA);
    pd = pd_traj(date, current_time, pd_raw, pd_noise, T, NoShunt, Pred_length, baseMVA)
    # ################################################################################
    # # for ploting demands
    # ################################################################################
    # # println(sum(pd.traj, dims=2))
    # # println(sum(pd.traj))
    # # println(size(pd.traj))
    # # println(size(pd.sigma))
    # # println(size(pd.ct))
    # plot(1:288, reshape(sum(pd.traj, dims=1), 288,1), label="traj")
    # plot(1:288, reshape(sum(pd.qd_traj, dims=1), 288,1), label="qd traj")
    # plot(1:288, reshape(sum(pd.sigma, dims=1), 288,1), label="sigma")

    # # println(maximum(sum(pd.traj, dims=1)))
    # # plot(1:288, reshape(pd.traj[1,:], 288,1), label="traj")
    # # plot!(1:288, reshape(pd.qd_traj[1,:], 288,1), label="qd traj")
    # ################################################################################
    # pg = pg_traj_pu_det(current_time, pg_raw, p_rate, T, NoShunt, baseMVA);
    pg = pg_traj(date, current_time, pg_raw, pg_noise, solar_error_max, p_rate, T, Pred_length,
    NoShunt, baseMVA);
    # # ################################################################################
    # # # for ploting solars
    # # ################################################################################
    # # # println(size(pg.mu))
    # # # println(size(pg.sigma))
    # # # println(size(pg.mu_ct))
    # # # println(pg.sg_max)
    # # # println(sum(pg.mu)/sum(pd.traj))
    # plot(1:288, reshape(sum(pg.mu, dims=1), 288,1), label="traj")
    # # ################################################################################
    #
    NoShunt = length(shunt_struct.P)
    NoBus = length(bus_struct.baseKV)

    Pd_bus = zeros(NoBus,1)
    Qd_bus = zeros(NoBus,1)
    for bus=1:NoBus
        bus_shunt_list = findall(id->id==bus, shunt_struct.find_bus[:,1]);
        if ~isempty(bus_shunt_list)
            Pd_bus[bus,1] = sum(
                pd.traj[Int(shunt),1]-pg.mu[Int(shunt),1]
                for shunt in bus_shunt_list);

            Qd_bus[bus,1] = sum(pd.qd_traj[Int(shunt),1]
                for shunt in bus_shunt_list);
        end
    end
    reference_points = reference_point(sum(Pd_bus), Pd_bus, Qd_bus);
    #
    #
    #
    obj = GML_Sys_Ava_NYISO(T, pd, ancillary_type, B_cap, icdf, shunt_struct,
        bus_struct, branch_struct, gen_struct, baseMVA);

    ############ write B and P_rsrv_feedback

    val_opt = optimal_NYISO(SN, current_time, obj, ancillary_type, baseMVA,
        feedback, pd, pg, price, shunt_struct, bus_struct, branch_struct, gen_struct,
        reference_points);

    B_feedback_out=zeros(NoShunt,1)
    for shunt = 1:NoShunt
        B_feedback_out[shunt, 1] = val_opt.B[shunt,1] -
            floor(val_opt.R[shunt,1]/12*baseMVA*1000)/1000;
    end
    write_B_out(current_time, B_feedback_out)
        # println(val_opt.P_rsrv)
    if ancillary_type == "10min" || ancillary_type == "30min"
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
        write_RSRV_out(P_rsrv_feedback_temp)
    end


    Pd_bus_real = zeros(NoBus,1)
    Qd_bus_real = zeros(NoBus,1)
    for bus=1:NoBus
        bus_shunt_list = findall(id->id==bus, shunt_struct.find_bus[:,1]);
        if ~isempty(bus_shunt_list)
            Pd_bus_real[bus,1] = sum(
                pd.ct[Int(shunt),1]-pg.mu_ct[Int(shunt),1]+val_opt.P_cul[Int(shunt),1]
                -val_opt.R[Int(shunt),1]
                for shunt in bus_shunt_list);

            Qd_bus_real[bus,1] = sum(pd.qd_ct[Int(shunt),1]-val_opt.Qg[Int(shunt),1]
                for shunt in bus_shunt_list);
        end
    end
    val_rp = reference_point(sum(Pd_bus), Pd_bus_real, Qd_bus_real);
    #
    mkpath("./result")
    write_output_out_2(val_opt, 0, val_rp, price, Pd_bus_real,
        string("./result/Time", current_time, ".csv"));
    write_v_bus_out(current_time, val_opt);
    return val_opt
end



function GML_Sys_Ava_NYISO(T, pd, ancillary_type, B_cap, icdf, shunt_struct,
    bus_struct, branch_struct, gen_struct, baseMVA)
    println("===== GML - Boundaries Buildup");
    ###############################################################################
    NoShunt = length(shunt_struct.P)
    NoBus = length(bus_struct.baseKV)
    NoBr = length(branch_struct.r)
    Qf_max = zeros(NoShunt,T)
    Qf_min = zeros(NoShunt,T)
    for shunt = 1:NoShunt
        # println(bus)
        if shunt_struct.type[shunt] == 2
            Qf_max[shunt, :]=shunt_struct.Q_max[shunt, 1]*ones(1,T)/baseMVA;
            # println(Qf_max[shunt, :])
        end
    end


    # minimum solar
    Pg_min = zeros(NoShunt, T);
    other_shunt = vcat(1, 3, 6:NoShunt);
    B_max_temp = reshape(shunt_struct.frac_P*B_cap, NoShunt, 1)*ones(1, T);
    B_ohter_shunt_sum = sum(B_max_temp[shunt,1] for shunt in other_shunt);
    B_max = B_max_temp*(B_cap-1.558)/B_ohter_shunt_sum;
    # println(sum(B_max)/288)
    B_min = zeros(NoShunt, T);
    R_rate = (375-7.382)/(75-1.558);
    R_max = R_rate*B_max/baseMVA;
    R_max[2,:] = 1.450/baseMVA*ones(1, T);
    R_max[4,:] = 2.265/baseMVA*ones(1, T);
    R_max[5,:] = 3.667/baseMVA*ones(1, T);
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
    V_min = 0.95;
    V_max = 1.05;
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
    feedback, pd, pg, price, shunt_struct, bus_struct, branch_struct, gen_struct,
    reference_points);
    # Q_gamma=0.01;
    Q_bar=0.01;
    P_Percent = [0.181935169, 0.002855468, 0.001789426, 0.004736269, 0.030496394,
        0.002840238, 0.000639625, 0.386858757, 0.387848653];
    # Q_Percent = [0.475276323, 0.034904014, -0.017452007, -0.029086678,
    #     -0.034904014, -0.031297266, -0.053868528, 0.357998837, 0.298429319];
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
    Qd = pd.qd_traj;

    Pd_agg = sum(Pd, dims=1);
    Qd_agg = sum(Qd, dims=1);

    B_feedback = feedback.B_feedback;
    P_rsrv_feedback = feedback.P_rsrv_feedback;

    NoShunt = length(shunt_struct.P)
    NoBus = length(bus_struct.baseKV)
    NoBr = length(branch_struct.fbus)
    NoGen = length(gen_struct.Pmax)
    bus_switch_shunt=shunt_struct.bus_switch_shunt;
    # NoGen=length(gen_data.id)

    # m = Model(with_optimizer(Mosek.Optimizer, QUIET=true,
    # MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-3,
    # MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-3))
    m = Model(Mosek.Optimizer)
    set_optimizer_attribute(m, "QUIET", true)
    # define the real-time variables

    # Peak-shaving
    @variable(m, 100 >=peak_withdraw[1, 1] >= 0)

    # Shunt level
    @variable(m, pg.mu[shunt,1] >= Pg_rt[shunt=1:NoShunt, 1] >= 0)
    @variable(m, pg.sg_max[shunt, 1] >= Qg_rt[shunt=1:NoShunt, 1] >= -pg.sg_max[shunt, 1])
    # @variable(m, temp[1:NoShunt, 1])
    @variable(m, B_max[shunt,1] >= B_rt[shunt = 1:NoShunt, 1] >= 0)
    @variable(m, R_max[shunt,1]>= R_rt[shunt = 1:NoShunt, 1] >= -R_max[shunt, 1])

    @variable(m, Pd[shunt,1]+R_max[shunt,1]>= P_shunt_rt[shunt = 1:NoShunt, 1]
        >= Pd[shunt,1]-pg.mu[shunt,1]-R_max[shunt,1])

    # @variable(m, P_shunt_rt[1:NoShunt, 1])
    # @variable(m, Qd[shunt,1]+pg.sg_max[shunt, 1] >= Q_shunt_rt[shunt = 1:NoShunt, 1]
    #    >=Qd[shunt,1] - pg.sg_max[shunt, 1])
    @variable(m, 1000>= Q_shunt_rt[1:NoShunt, 1] >=-1000)

    # bus level
    @variable(m, 1000 >= P_bus_rt[1:NoBus, 1]  >= -1000)
    @variable(m, 1000 >= Q_bus_rt[1:NoBus, 1] >= -1000)

    @variable(m, v_rt[1:NoBus, 1])

    # Branch level
    @variable(m, 1000 >= P_br_rt[1:NoBr, 1] >= -1000)
    @variable(m, 1000 >= Q_br_rt[1:NoBr, 1] >= -1000)
    @variable(m, l_rt[1:NoBr, 1])

    # gen level
    @variable(m, 1000 >= P_gen_rt[1:NoGen, 1]  >= -1000)
    @variable(m, 1000 >= Q_gen_rt[1:NoGen, 1]  >= -1000)


    # println(" ---- Real Time Constraint Buildup ")
    for shunt=1:NoShunt
        # box constraints on Solar (Pg), Battery discharge (R)
        # @constraint(m, 0<=Pg_rt[shunt,1]);
        # @constraint(m, R_min[shunt,1]<= R_rt[shunt,1]);
        # @constraint(m, R_rt[shunt,1]<= R_max[shunt,1]);
        # Initial battery SOC
        if B_feedback[shunt,1]>= B_max[shunt,1]
            B_feedback[shunt,1]= B_max[shunt,1];
        elseif B_feedback[shunt,1]<= B_min[shunt,1]
            B_feedback[shunt,1]= B_min[shunt,1];
        end
        @constraint(m, B_rt[shunt,1]==B_feedback[shunt,1]);


        if shunt_struct.type[shunt, 1]==1
            # solar
            # @constraint(m, Pg_rt[shunt,1]<=
            #     positive_scalar(icdf*sqrt(pd.sigma[shunt,1]+pg.sigma[shunt,1])
            #     +pg.mu[shunt,1]));
            # @constraint(m, Qg_rt[shunt,1] <= pg.sg_max[shunt, 1])
            # @constraint(m, Qg_rt[shunt,1] >= -pg.sg_max[shunt, 1])
            @constraint(m, [pg.sg_max[shunt, 1], Pg_rt[shunt,1], Qg_rt[shunt,1]]
             in SecondOrderCone());


            # SOC constrains on real and reactive power on bus
            # maxS_rt = sqrt(1+Q_gamma^2)*(Pd[shunt,1]+R_max[shunt,1]-
            # positive_scalar(
            # icdf*sqrt(pd.sigma[shunt,1]+pg.sigma[shunt,1])));
            # if maxS_rt>=0
            #     @constraint(m,
            #     [maxS_rt,
            #     P_shunt_rt[shunt,1], Q_shunt_rt[shunt,1]] in SecondOrderCone())
            # else
            #     @constraint(m,
            #     [-maxS_rt,
            #     P_shunt_rt[shunt,1], Q_shunt_rt[shunt,1]] in SecondOrderCone())
            # end

            # real power on bus is demand minus solar and discharge
            # AKA power injection

            @constraint(m, P_shunt_rt[shunt,1]==
                Pd[shunt,1]
                -(Pg_rt[shunt,1]+R_rt[shunt,1]));
            @constraint(m, Q_shunt_rt[shunt,1]==
                Qd[shunt,1]-Qg_rt[shunt,1]);
            # @constraint(m, Q_shunt_rt[shunt,1]==
            #     Qd[shunt,1]);
        elseif shunt_struct.type[shunt, 1]==2
            @constraint(m, Pg_rt[shunt,1]==0)
            @constraint(m, Qg_rt[shunt,1]==0)
            @constraint(m, P_shunt_rt[shunt,1]==0);
            if Qd[shunt,1] >= Q_bar;
                @constraint(m, Q_shunt_rt[shunt,1]<=0);
                @constraint(m, Q_shunt_rt[shunt,1]>=-Qf_max[shunt,1]);
            else
                @constraint(m, Q_shunt_rt[shunt,1]==0);
            end
        end
    end



    for bus=1:NoBus
        # id of shunts that belong to that bus
        bus_shunt_list = findall(id->id==bus, shunt_struct.find_bus[:,1]);
        # a list of branchs that FROM the bus of interest
        sub_branch_list = findall(one->one==1, C_ind[bus,:]) # Out
        # a list of branchs that POINT TO the bus of interest
        add_branch_list = findall(minusone->minusone==-1, C_ind[bus,:]) #In

        if isempty(bus_shunt_list)
            @constraint(m, P_bus_rt[bus,1] == 0)
            @constraint(m, Q_bus_rt[bus,1] == 0)
        elseif ~isempty(bus_shunt_list)
            @constraint(m, P_bus_rt[bus,1] ==
                sum(P_shunt_rt[Int(shunt),1] for shunt in bus_shunt_list))
            @constraint(m, Q_bus_rt[bus,1] ==
                sum(Q_shunt_rt[Int(shunt),1] for shunt in bus_shunt_list))
        end

        ########### voltage
        # if bus == 70
        #     @constraint(m,
        #         1.03569^2 == v_rt[bus,1]);
        # else
        #     @constraint(m, V_min^2<= v_rt[bus,1]);
        #     @constraint(m, v_rt[bus,1]<= V_max^2);
        # end
        if bus_switch_shunt[bus] == 0
            @constraint(m, V_min^2<= v_rt[bus,1]);
            @constraint(m, v_rt[bus,1]<= V_max^2);
        else
            switch_shunt_id = Int(bus_switch_shunt[bus])
            @constraint(m,
                shunt_struct.Vsp[switch_shunt_id]^2 == v_rt[bus,1]);
        end

        if bus_struct.type[bus]==1 # non-generator bus

            # Power Balance Equations
            # Power injection = Power Flow In - Power Flow Out
            if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    sum(P_br_rt[branch,1] for branch in add_branch_list)
                    -sum(P_br_rt[branch,1] for branch in sub_branch_list)
                    );
                @constraint(m, Q_bus_rt[bus,1]==
                    sum(Q_br_rt[branch,1] for branch in add_branch_list)
                    -sum(x[branch]*l_rt[branch,1] for branch in add_branch_list)
                    -sum(Q_br_rt[branch,1] for branch in sub_branch_list));
            elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    sum(P_br_rt[branch,1] for branch in add_branch_list)
                    );
                @constraint(m, Q_bus_rt[bus,1]==
                    sum(Q_br_rt[branch,1] for branch in add_branch_list)
                    -sum(x[branch]*l_rt[branch,1] for branch in add_branch_list)
                    );
            elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    -sum(P_br_rt[branch,1] for branch in sub_branch_list)
                    );
                @constraint(m, Q_bus_rt[bus,1]==
                    -sum(Q_br_rt[branch,1] for branch in sub_branch_list));
            end
        else

            # identify the id of the generator that connects to the bus of interest
            gen_id = findall(bus_id ->bus_id == bus, gen_struct.bus[:,1]);
            # println(gen_id)
            # Power Balance Equations
            # Power injection = Power Flow In - Power Flow Out + Power Generatered
            if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    sum(P_br_rt[branch,1] for branch in add_branch_list)
                    -sum(P_br_rt[branch,1] for branch in sub_branch_list)
                    +P_gen_rt[gen_id[1], 1]);
                @constraint(m, Q_bus_rt[bus,1]==
                    sum(Q_br_rt[branch,1] for branch in add_branch_list)
                    -sum(x[branch]*l_rt[branch,1] for branch in add_branch_list)
                    -sum(Q_br_rt[branch,1] for branch in sub_branch_list)
                    +Q_gen_rt[gen_id[1], 1]);
            elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    sum(P_br_rt[branch,1] for branch in add_branch_list)
                    +P_gen_rt[gen_id[1], 1]);
                @constraint(m, Q_bus_rt[bus,1]==
                    sum(Q_br_rt[branch,1] for branch in add_branch_list)
                    -sum(x[branch]*l_rt[branch,1] for branch in add_branch_list)
                    +Q_gen_rt[gen_id[1], 1]);
            elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
                @constraint(m, P_bus_rt[bus,1]==
                    -sum(P_br_rt[branch,1] for branch in sub_branch_list)
                    +P_gen_rt[gen_id[1], 1]);
                @constraint(m, Q_bus_rt[bus,1]==
                    -sum(Q_br_rt[branch,1] for branch in sub_branch_list)
                    +Q_gen_rt[gen_id[1], 1]);
            end
        # elseif bus_struct.type[bus]==3
        #
        #     # identify the id of the generator that connects to the bus of interest
        #     gen_id = findall(bus_id ->bus_id == bus, gen_struct.bus[:,1]);
        #     # println(gen_id)
        #     # Power Balance Equations
        #     # Power injection = Power Flow In - Power Flow Out + Power Generatered
        #     if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
        #         @constraint(m, P_bus_rt[bus,1]==
        #             sum(P_br_rt[branch,1] for branch in add_branch_list)
        #             -sum(P_br_rt[branch,1] for branch in sub_branch_list)
        #             +P_gen_rt[gen_id[1], 1]+P_gen_rt[NoGen+1, 1]);
        #         @constraint(m, Q_bus_rt[bus,1]==
        #             sum(Q_br_rt[branch,1] for branch in add_branch_list)
        #             -sum(Q_br_rt[branch,1] for branch in sub_branch_list)
        #             +Q_gen_rt[gen_id[1], 1]+Q_gen_rt[NoGen+1, 1]);
        #     elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
        #         @constraint(m, P_bus_rt[bus,1]==
        #             sum(P_br_rt[branch,1] for branch in add_branch_list)
        #             +P_gen_rt[gen_id[1], 1]+P_gen_rt[NoGen+1, 1]);
        #         @constraint(m, Q_bus_rt[bus,1]==
        #             sum(Q_br_rt[branch,1] for branch in add_branch_list)
        #             +Q_gen_rt[gen_id[1], 1]+Q_gen_rt[NoGen+1, 1]);
        #     elseif isempty(add_branch_list) && ~isempty(sub_branch_list)
        #         @constraint(m, P_bus_rt[bus,1]==
        #             -sum(P_br_rt[branch,1] for branch in sub_branch_list)
        #             +P_gen_rt[gen_id[1], 1]+P_gen_rt[NoGen+1, 1]);
        #         @constraint(m, Q_bus_rt[bus,1]==
        #             -sum(Q_br_rt[branch,1] for branch in sub_branch_list)
        #             +Q_gen_rt[gen_id[1], 1]+Q_gen_rt[NoGen+1, 1]);
        #     end
        end
    end
    # voltage constraint for LinDistFlow
    for branch =1:NoBr
        fbus = branch_struct.fbus[branch];
        tbus = branch_struct.tbus[branch];

        @constraint(m, v_rt[fbus,1]-v_rt[tbus,1]==
            2*(r[branch]*P_br_rt[branch,1]+x[branch]*Q_br_rt[branch,1])
            -(r[branch]^2+x[branch]^2)*l_rt[branch,1]);

        # @constraint(m, 2*reference_points.p_star[branch,1]*P_br_rt[branch]
        # +2*reference_points.q_star[branch,1]*Q_br_rt[branch]
        # -reference_points.v_star[fbus,1]*l_rt[branch,1]
        # -reference_points.l_star*v_rt[fbus,1] == 0
        # );
        @constraint(m, 2*reference_points.p_star[branch,1]*P_br_rt[branch,1]
        +2*reference_points.q_star[branch,1]*Q_br_rt[branch,1]
        -reference_points.v_star[fbus,1]*l_rt[branch,1]
        -reference_points.l_star[branch,1]*v_rt[fbus,1]
        == 0
        );
    end

    # # box constraint on generator
    # # NoGen
    for Gen=2:NoGen
        # @constraint(m,
        #     P_gen_rt[Gen,1]<=gen_struct.Pmax[Gen]/baseMVA)
        # @constraint(m,
        #     P_gen_rt[Gen,1]>=gen_struct.Pmin[Gen]/baseMVA)
        @constraint(m,
            P_gen_rt[Gen,1]==Pd_agg[1,1]*P_Percent[Gen])
        # @constraint(m,
        #     Q_gen_rt[Gen,1]==Qd_agg[1,1]*Q_Percent[Gen])
        @constraint(m,
            Q_gen_rt[Gen,1]<=gen_struct.Qmax[Gen]/baseMVA)
        @constraint(m,
            Q_gen_rt[Gen,1]>=gen_struct.Qmin[Gen]/baseMVA)
    end
    # peak shaving
    @constraint(m,
            P_gen_rt[1,1]<=peak_withdraw[1, 1]);
    @constraint(m,
            1000>=peak_withdraw[1, 1]);


    #
    # # shunt
    @variable(m, pg.mu[shunt, t+1]>=
        Pg[1:SN, shunt = 1:NoShunt, t=1:T-1]>= 0); # the real power output
    @variable(m, pg.sg_max[shunt, 1]>= Qg[1:SN, shunt = 1:NoShunt, 1:T-1]
        >= -pg.sg_max[shunt, 1]); # the real power output

    @variable(m, B_max[shunt,t+1] >=B[1:SN, shunt=1:NoShunt, t=1:T-1] >= 0); # the storage
    @variable(m, R_max[shunt, t+1]>=
        R[1:SN, shunt=1:NoShunt, t=1:T-1] >= R_min[shunt, t+1]);# the charge/discharge rate

    @variable(m, Pd[shunt,t+1]+R_max[shunt, t+1] >=
        P_shunt[1:SN, shunt=1:NoShunt, t=1:T-1] >=
        Pd[shunt,t+1]-pg.mu[shunt, t+1]-R_max[shunt, t+1])
    @variable(m, 1000>=Q_shunt[1:SN, 1:NoShunt, 1:T-1]>=-1000)


    @variable(m, 1000>=P_bus[1:SN, 1:NoBus, 1:T-1]>=-1000)
    @variable(m, 1000>=Q_bus[1:SN, 1:NoBus, 1:T-1]>=-1000)
    @variable(m, v[1:SN, 1:NoBus, 1:T-1])

    # Branch level
    @variable(m, 1000>=P_br[1:SN,1:NoBr, 1:T-1]>=-1000)
    @variable(m, 1000>=Q_br[1:SN,1:NoBr, 1:T-1]>=-1000)
    @variable(m, l[1:SN,1:NoBr, 1:T-1])
    # Gen Level
    @variable(m, 1000>=P_gen[1:SN, 1:NoGen, 1:T-1]>=-1000)
    @variable(m, 1000>=Q_gen[1:SN, 1:NoGen, 1:T-1]>=-1000)


    # # P0
    # @variable(m, P0[1:SN, 1, 1:T-1])
    # @variable(m, Q0[1:SN, 1, 1:T-1])
    #
    # # for different scenario
    for scenario = 1:SN
        # for different time at prediction horizion
        for t=1:T-1
            # for different bus
            for shunt=1:NoShunt
                # box constraints on Solar (Pg), Battery discharge (R)
                # @constraint(m, 0<=Pg[scenario, shunt ,t]);
                # @constraint(m, R_min[shunt,t+1] <= R[scenario,shunt,t]);
                # @constraint(m, R[scenario,shunt,t]<= R_max[shunt,t+1]);
                # battery box constraint on SOC
                # @constraint(m, B_min[shunt,t+1] <= B[scenario,shunt,t]);
                # @constraint(m, B[scenario,shunt,t]<= B_max[shunt,t+1]);
                # battery SOC dynamics
                if t==1
                    @constraint(m, B[scenario,shunt,1] == B_rt[shunt,1]
                        -delta_t*(R_rt[shunt,1]*baseMVA))
                else
                    @constraint(m, B[scenario,shunt,t] ==
                        B[scenario,shunt,t-1] - R[scenario,shunt,t-1]*baseMVA*delta_t)
                end
                if shunt_struct.type[shunt, 1]==1
                    # @constraint(m, Pg[scenario, shunt,t]<=
                    #     positive_scalar(
                    #     icdf*sqrt(pd.sigma[shunt,t+1]+pg.sigma[shunt,t+1])+pg.mu[shunt,t+1])
                    #     );
                    # @constraint(m, Qg[scenario, shunt,t]<=pg.sg_max[shunt, 1])
                    # @constraint(m, Qg[scenario, shunt,t]>=-pg.sg_max[shunt, 1])
                    @constraint(m,
                        [ pg.sg_max[shunt, 1], Pg[scenario, shunt,t], Qg[scenario, shunt,t] ]
                        in SecondOrderCone());
                    # @constraint(m,
                    #     Qg[scenario, shunt, t]<=pg.sg_max[shunt, 1]);
                    # real power on bus is demand minus solar and discharge
                    # AKA power injection
                    @constraint(m,  P_shunt[scenario,shunt,t]==
                        Pd[shunt,t+1]
                        -Pg[scenario,shunt,t]-R[scenario, shunt,t]);
                    @constraint(m,  Q_shunt[scenario,shunt,t]==
                        Qd[shunt,t+1]-Qg[scenario,shunt,t]);

                    # SOC constrains on real and reactive power on bus
                    # maxS = sqrt(1+Q_gamma^2)*(Pd[shunt,t+1]-positive_scalar(
                    # icdf*sqrt(pd.sigma[shunt,t+1]+pg.sigma[shunt,t+1])+pg.mu[shunt,t+1])
                    # -R_max[shunt,t+1]);
                    # if maxS >=0
                    #     @constraint(m,
                    #     [maxS, P_shunt[scenario,shunt,t], Q_shunt[scenario,shunt,t]]
                    #      in SecondOrderCone())
                    # else
                    #     @constraint(m,
                    #     [-maxS, P_shunt[scenario,bus,t], Q_shunt[scenario,bus,t]]
                    #      in SecondOrderCone())
                    # end
                elseif shunt_struct.type[shunt, 1]==2
                    @constraint(m, Pg[scenario, shunt, t]==0)
                    @constraint(m, Qg[scenario,shunt,t]==0)
                    @constraint(m, P_shunt[scenario,shunt,t]==0)

                    if Qd[shunt, t+1] >= Q_bar;
                        @constraint(m, Q_shunt[scenario,shunt,t]<=0)
                        @constraint(m, Q_shunt[scenario,shunt,t]>=-Qf_max[shunt,1])
                    else
                        @constraint(m, Q_shunt[scenario,shunt,t]==0)
                    end

                end
            end
            for bus=1:NoBus
                # id of shunts that belong to that bus
                bus_shunt_list = findall(id->id==bus, shunt_struct.find_bus[:,1]);
                # a list of branchs that FROM the bus of interest
                sub_branch_list = findall(one->one==1, C_ind[bus,:]) # out
                # a list of branchs that POINT TO the bus of interest
                add_branch_list = findall(minusone->minusone==-1, C_ind[bus,:]) # In

                if isempty(bus_shunt_list)
                    @constraint(m, P_bus[scenario, bus, t] == 0)
                    @constraint(m, Q_bus[scenario, bus, t] == 0)
                elseif ~isempty(bus_shunt_list)
                    @constraint(m, P_bus[scenario, bus, t] ==
                        sum(P_shunt[scenario,Int(shunt),t] for shunt in bus_shunt_list))
                    @constraint(m, Q_bus[scenario, bus, t] ==
                        sum(Q_shunt[scenario,Int(shunt),t] for shunt in bus_shunt_list))
                end

                ######## voltage
                if bus_switch_shunt[bus] == 0
                    # box constraint on bus
                    @constraint(m, V_min^2<= v[scenario,bus,t]);
                    @constraint(m, v[scenario,bus,t]<= V_max^2);
                else
                    switch_shunt_id = Int(bus_switch_shunt[bus])
                    @constraint(m,
                        shunt_struct.Vsp[switch_shunt_id]^2== v[scenario,bus,t]);
                end
                # if bus == 70
                #     @constraint(m,
                #         1.03569^2 == v[scenario,bus,t]);
                # else
                #     @constraint(m, V_min^2<= v[scenario,bus,t]);
                #     @constraint(m, v[scenario,bus,t]<= V_max^2);
                # end

                if bus_struct.type[bus] == 1 # non-generator bus

                    # Power Balance Equations
                    # Power injection = Power Flow In - Power Flow Out
                    if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            sum(P_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
                            );
                        @constraint(m, Q_bus[scenario,bus,t]==
                            sum(Q_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(x[branch]*l[scenario,branch,t] for branch in add_branch_list)
                            -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
                            );
                    elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            sum(P_br[scenario,branch,t] for branch in add_branch_list)
                            );
                        @constraint(m, Q_bus[scenario,bus,t]==
                            sum(Q_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(x[branch]*l[scenario,branch,t] for branch in add_branch_list)
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

                    gen_id = findall(bus_id ->bus_id == bus, gen_struct.bus[:,1]);

                    if ~isempty(add_branch_list) && ~isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            sum(P_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(P_br[scenario,branch,t] for branch in sub_branch_list)
                            +P_gen[scenario,gen_id[1],t]);
                        @constraint(m, Q_bus[scenario,bus,t]==
                            sum(Q_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(x[branch]*l[scenario,branch,t] for branch in add_branch_list)
                            -sum(Q_br[scenario,branch,t] for branch in sub_branch_list)
                            +Q_gen[scenario,gen_id[1],t]);
                    elseif ~isempty(add_branch_list) && isempty(sub_branch_list)
                        @constraint(m, P_bus[scenario,bus,t]==
                            sum(P_br[scenario,branch,t] for branch in add_branch_list)
                            +P_gen[scenario,gen_id[1],t]);
                        @constraint(m, Q_bus[scenario,bus,t]==
                            sum(Q_br[scenario,branch,t] for branch in add_branch_list)
                            -sum(x[branch]*l[scenario,branch,t] for branch in add_branch_list)
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
                    +x[branch]*Q_br[scenario,branch,t])
                    -(r[branch]^2+x[branch]^2)*l[scenario,branch,t]);

                @constraint(m,
                    2*reference_points.p_star[branch,1]*P_br[scenario,branch,t]
                    +2*reference_points.q_star[branch,1]*Q_br[scenario,branch,t]
                    -reference_points.v_star[fbus,1]*l[scenario,branch,t]
                    -reference_points.l_star[branch,1]*v[scenario,fbus,t]
                    == 0
                    );
            end
            for Gen=2:NoGen
                @constraint(m,
                    P_gen[scenario,Gen,t]==Pd_agg[1,t+1]*P_Percent[Gen])
                # @constraint(m,
                #     Q_gen[scenario,Gen,t]==Qd_agg[1,t+1]*Q_Percent[Gen])
                # @constraint(m,
                #     P_gen[scenario,Gen,t]<=gen_struct.Pmax[Gen]/baseMVA)
                # @constraint(m,
                #     P_gen[scenario,Gen,t]>=gen_struct.Pmin[Gen]/baseMVA)
                @constraint(m,
                    Q_gen[scenario,Gen,t]<=gen_struct.Qmax[Gen]/baseMVA)
                @constraint(m,
                    Q_gen[scenario,Gen,t]>=gen_struct.Qmin[Gen]/baseMVA)
            end
            @constraint(m,
                P_gen[scenario,1,t]<=peak_withdraw[1, 1]);
            # @constraint(m,
            #     1000>=peak_withdraw[1, 1]);
        end
    end
    #
    #
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
            P_rsrv,price,pg,pd,beta,SN,obj, gen_struct, bus_struct,
            branch_struct, peak_withdraw,baseMVA))
    else
        @objective(m, Min,
            fn_cost_RHC_rt(delta_t, P_gen, P_gen_rt, Pg_rt,Pg, price,
                pg, pd, beta, SN, obj, gen_struct, bus_struct,
                branch_struct, peak_withdraw, baseMVA))
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
    println(string("    ----Optimal Cost of the whole traj: ", cost_o))

    println(string("    ----current demand: ", baseMVA*sum(Pd[:,1])))


    probability = price.probability[1:SN];
    scenario_of_interest = findall(ind->ind==maximum(probability), probability)[1]
    # println(scenario_of_interest);

    ###############
    R_rt_o=JuMP.value.(R_rt)
    R_o=JuMP.value.(R)
    R_sum_ct = sum(R_rt_o[:, 1])
    println(string("    ----current battery disharge: ", R_sum_ct))
    R_traj = hcat(R_rt_o, R_o[scenario_of_interest,:,:])
    ###############
    P_br_rt_o=JuMP.value.(P_br_rt)
    P_br_o=JuMP.value.(P_br)
    P_br_traj = hcat(P_br_rt_o, P_br_o[scenario_of_interest,:,:])

    ###############
    Q_br_rt_o=JuMP.value.(Q_br_rt)
    Q_br_o=JuMP.value.(Q_br)
    Q_br_traj = hcat(Q_br_rt_o, Q_br_o[scenario_of_interest,:,:])
    ###############
    v_rt_o=JuMP.value.(v_rt)
    v_o=JuMP.value.(v)
    v_traj = hcat(v_rt_o, v_o[scenario_of_interest,:,:])
    ################
    B_rt_o=JuMP.value.(B_rt)
    B_o=JuMP.value.(B)
    B_traj = hcat(B_rt_o, B_o[scenario_of_interest,:,:])
    ###############
    Pg_rt_o=JuMP.value.(Pg_rt)
    Pg_o=JuMP.value.(Pg)
    Pg_traj = hcat(Pg_rt_o, Pg_o[scenario_of_interest,:,:])
    #

    Qg_rt_o=JuMP.value.(Qg_rt)
    Qg_o=JuMP.value.(Qg)
    Qg_traj = hcat(Qg_rt_o, Qg_o[scenario_of_interest,:,:])

    ############
    P_shunt_rt_o=JuMP.value.(P_shunt_rt)
    P_shunt_o=JuMP.value.(P_shunt)
    P_shunt_traj = hcat(P_shunt_rt_o, P_shunt_o[scenario_of_interest,:,:])
    # println(string("    ----current real on shunt: ",
    #     sum(P_shunt_traj)))

    ############
    Q_shunt_rt_o=JuMP.value.(Q_shunt_rt)
    Q_shunt_o=JuMP.value.(Q_shunt)
    Q_shunt_traj = hcat(Q_shunt_rt_o, Q_shunt_o[scenario_of_interest,:,:])
    # println(string("    ----current reactive on shunt: ",
    #     sum(Q_shunt_traj)))
    ############
    P_bus_rt_o=JuMP.value.(P_bus_rt)
    P_bus_o=JuMP.value.(P_bus)
    P_bus_traj = hcat(P_bus_rt_o, P_bus_o[scenario_of_interest,:,:])

    ############
    Q_bus_rt_o=JuMP.value.(Q_bus_rt)
    Q_bus_o=JuMP.value.(Q_bus)
    Q_bus_traj = hcat(Q_bus_rt_o, Q_bus_o[scenario_of_interest,:,:])
    println(string("    ----current reactive on bus: ",
        sum(Q_bus_traj)))

    ############
    P_gen_rt_o=JuMP.value.(P_gen_rt)
    P_gen_o=JuMP.value.(P_gen)
    P_gen_traj = hcat(P_gen_rt_o, P_gen_o[scenario_of_interest,:,:])
    P_gen_sum_traj = hcat(sum(P_gen_rt_o), sum(P_gen_o[scenario_of_interest,:,:], dims=1))
    # println(size(P_gen_sum_traj))

    P_gen_sum_ct = sum(P_gen_rt_o[gen, 1]*baseMVA for gen in 1:NoGen)
    println(string("    ----current generation: ", P_gen_sum_ct))


    ############
    Q_gen_rt_o = JuMP.value.(Q_gen_rt)
    Q_gen_o = JuMP.value.(Q_gen)
    Q_gen_traj = hcat(Q_gen_rt_o, Q_gen_o[scenario_of_interest,:,:])
    Q_gen_sum_traj = hcat(sum(Q_gen_rt_o), sum(Q_gen_o[scenario_of_interest,:,:], dims=1))
    Q_gen_sum_ct = sum(Q_gen_rt_o[gen, 1]*baseMVA for gen in 1:NoGen)
    println(string("    ----current reactive generation: ",
        sum(Q_gen_sum_ct)))

    ############
    peak_withdraw_o = JuMP.value.(peak_withdraw);
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
        P_rsrv_total = hcat(P_rsrv_rt_o[1,1], reshape(P_rsrv_s[scenario_of_interest,:], 1, 287));
        B_rsrv_rt_o=JuMP.value(B_rsrv_rt);
        B_rsrv_s=JuMP.value.(B_rsrv);
        B_rsrv_total = hcat(B_rsrv_rt_o[1,1], reshape(B_rsrv_s[scenario_of_interest,:], 1, 287));
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
    P_cul = zeros(length(pg.mu[:,1]),1)
    for shunt = 1: length(pg.mu[:,1])
        P_cul[shunt,1] = pg.mu[shunt,1]-Pg_rt_o[shunt,1]
    end
    println(string("    ----curtailment solar: ", sum(P_cul))
    alpha_1 = price.alpha_scenario[1,:]
    lambda_1 = price.lambda_scenario[1,:]

    println(string("    ----Optimzal cost at this instance: ", Cost_real))
    pg_upper = hcat(
    icdf*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])+pg.mu[:,:])

    ########################
    # feeders of interest
    FOL_pred_length = 12;

    R_FOL = zeros(3, FOL_pred_length)
    Pg_FOL = zeros(3, FOL_pred_length)
    P_FOL = zeros(3, FOL_pred_length)
    Q_FOL = zeros(3, FOL_pred_length)
    v_FOL_sqr = zeros(3, FOL_pred_length)
    v_FOL = zeros(3, FOL_pred_length)
    Feeder_shunt_number = [4, 2, 5]
    Feeder_bus_number = [3, 2, 3]
    scenario_poll = 1:SN

    for feeder = 1:3
        R_FOL[feeder, :] = hcat(R_rt_o[Feeder_shunt_number[feeder], 1],
            reshape(sum(probability[scenario]
            *R_o[scenario,Feeder_shunt_number[feeder],1:FOL_pred_length-1] for scenario in scenario_poll)
            , 1, FOL_pred_length-1));
        Pg_FOL[feeder, :] = hcat(Pg_rt_o[Feeder_shunt_number[feeder], 1],
            reshape(sum(probability[scenario]
            *Pg_o[scenario,Feeder_shunt_number[feeder],1:FOL_pred_length-1] for scenario in scenario_poll)
            , 1, FOL_pred_length-1));
        P_FOL[feeder, :] = hcat(P_shunt_rt_o[Feeder_shunt_number[feeder], 1],
            reshape(sum(probability[scenario]
            *P_shunt_o[scenario,Feeder_shunt_number[feeder],1:FOL_pred_length-1] for scenario in scenario_poll)
            , 1, FOL_pred_length-1));
        Q_FOL[feeder, :] = hcat(Q_shunt_rt_o[Feeder_shunt_number[feeder], 1],
            reshape(sum(probability[scenario]
            *Q_shunt_o[scenario,Feeder_shunt_number[feeder],1:FOL_pred_length-1] for scenario in scenario_poll)
            , 1, FOL_pred_length-1));
        v_FOL_sqr[feeder, :] = hcat(v_rt_o[Feeder_bus_number[feeder], 1],
            reshape(sum(probability[scenario]
            *v_o[scenario,Feeder_bus_number[feeder],1:FOL_pred_length-1] for scenario in scenario_poll)
            , 1, FOL_pred_length-1));
    end
    v_FOL=sqrt.(v_FOL_sqr)


    # R_FOL = vcat(
    #     reshape(sum(R_traj[2,1:FOL_pred_length]),1,FOL_pred_length),
    #     reshape(R_traj[4,1:FOL_pred_length],1,FOL_pred_length),
    #     reshape(R_traj[5,1:FOL_pred_length],1,FOL_pred_length));
    # Pg_FOL = vcat(
    #         reshape(Pg_traj[2,1:FOL_pred_length],1,FOL_pred_length),
    #         reshape(Pg_traj[4,1:FOL_pred_length],1,FOL_pred_length),
    #         reshape(Pg_traj[5,1:FOL_pred_length],1,FOL_pred_length));
    # Q_FOL = vcat(
    #     reshape(Q_shunt_traj[2,1:FOL_pred_length],1,FOL_pred_length),
    #     reshape(Q_shunt_traj[4,1:FOL_pred_length],1,FOL_pred_length),
    #     reshape(Q_shunt_traj[5,1:FOL_pred_length],1,FOL_pred_length));
    # P_FOL = vcat(
    #     reshape(P_shunt_traj[2,1:FOL_pred_length],1,FOL_pred_length),
    #     reshape(P_shunt_traj[4,1:FOL_pred_length],1,FOL_pred_length),
    #     reshape(P_shunt_traj[5,1:FOL_pred_length],1,FOL_pred_length));
    # v_FOL = vcat(
    #     reshape(sqrt.(v_traj[1,1:FOL_pred_length]),1,FOL_pred_length),
    #     reshape(sqrt.(v_traj[2,1:FOL_pred_length]),1,FOL_pred_length),
    #     reshape(sqrt.(v_traj[2,1:FOL_pred_length]),1,FOL_pred_length));

    ########################

    val_opt = (R=(R_rt_o), B=(B_rt_o), Pg=(Pg_rt_o), Qg=(Qg_rt_o),
        P_gen=(P_gen_rt_o), Q_gen=(Q_gen_rt_o),
        P_rsrv=(P_rsrv_rt_o), B_rsrv=(B_rsrv_rt_o),
        P_bus=(P_bus_rt_o), Q_bus=(Q_bus_rt_o),
        v=(v_rt_o), P_br=(P_br_rt_o), Q_br=(Q_br_rt_o),
        Cost_real=(Cost_real), time_solve=(time_solve),
        P_rsrv_total=(P_rsrv_total), B_rsrv_total=(B_rsrv_total),
        alpha_1=(alpha_1),
        R_traj = (R_traj), B_traj=(B_traj),
        Pg_traj=(Pg_traj), Qg_traj=(Qg_traj),
        P_gen_traj=(P_gen_traj),
        Q_gen_traj=(Q_gen_traj),
        P0_traj = (P_gen_traj), P_bus_traj=(P_bus_traj), Q_bus_traj=(Q_bus_traj),
        pg_upper=(pg_upper),
        lambda_1=(lambda_1), pd=(pd), pg=(pg), v_traj=(v_traj), P_br_traj=(P_br_traj),
        Q_br_traj=(Q_br_traj), terminate_s = (terminate_s), P_cul = (P_cul),
        R_FOL = (R_FOL), Pg_FOL = (Pg_FOL), Q_FOL = (Q_FOL),
        P_FOL = (P_FOL), v_FOL = (v_FOL), peak_withdraw = (peak_withdraw_o),
        cost_agg = (cost_o-100*10000*peak_withdraw_o[1,1]));


    return val_opt
end


function fn_cost_RHC_rt(delta_t, P_gen, P_gen_rt, Pg_rt, Pg, price,
    pg, pd, beta, SN, obj, gen_struct, bus_struct,
    branch_struct, peak_withdraw, baseMVA)

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

    ################### Peak Shaving###################

    # Cost_peakshaving = sum(peak_withdraw[1,1])*baseMVA*10000;
    Cost_peakshaving = peak_withdraw[1,1]*baseMVA*10000;

    Final_cost = (Cost_P_gen_sum_ct+Cost_P_gen_sum_scenario
        +Cost_Pg_ct_diff+Cost_Pg_scenario_diff
        +Cost_peakshaving)

    return Final_cost
end

function fn_cost_RHC_anc(delta_t, P_gen, P_gen_rt, Pg_rt,Pg, P_rsrv_rt,
    P_rsrv,price,pg,pd,beta,SN,obj, gen_struct, bus_struct,
    branch_struct, peak_withdraw, baseMVA)

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

    ################### Peak Shaving###################

    Cost_peakshaving = peak_withdraw[1,1]*baseMVA*10000;

    Final_cost = (Cost_P_gen_sum_ct+Cost_P_gen_sum_scenario
        +Cost_Pg_ct_diff+Cost_Pg_scenario_diff
        -Revenue_P_rsrv_ct-Revenue_P_rsrv_scenario
        +Cost_peakshaving)
    # Final_cost = (Cost_P_gen_sum_ct)
    return Final_cost
end










function write_branch_real_output(t,val_opt)
    output = hcat(val_opt.P_br_traj)
    filename = string("t=",t,"/branch_flow_real_power_output.csv")
    CSV.write(filename, DataFrame(output));
    println("    ---- Finish branch flow real files! ")
end

function write_branch_reactive_output(t,val_opt)
    output = hcat(val_opt.Q_br_traj)
    filename = string("t=",t,"/branch_flow_reactive_power_output.csv")
    CSV.write(filename, DataFrame(output));
    println("    ---- Finish branch flow reactive files! ")
end

function write_bus_real_output(t,val_opt)
    output = hcat(val_opt.P_bus_traj)
    filename = string("t=",t,"/bus_real_output.csv")
    CSV.write(filename, DataFrame(output));
    println("    ---- Finish bus real files! ")
end

function write_bus_reactive_output(t,val_opt)
    output = hcat(val_opt.Q_bus_traj)
    filename = string("t=",t,"/bus_reactive_output.csv")
    CSV.write(filename, DataFrame(output));
    # CSV.write("bus_reactive_output.csv", DataFrame(output));
    println("    ---- Finish bus reactive files! ")
end

function write_bus_voltage_output(t,val_opt)
    output = hcat(sqrt.(val_opt.v_traj))
    filename = string("t=",t,"/bus_voltage_output.csv")
    CSV.write(filename, DataFrame(output));
    # CSV.write("bus_voltage_output.csv", DataFrame(output));
    println("    ---- Finish bus voltage files! ")
end

function write_generator_real_output(t,val_opt)
    output = hcat(val_opt.P_gen_traj)
    filename = string("t=",t,"/gen_real_output.csv")
    CSV.write(filename, DataFrame(output));
    # CSV.write("gen_real_output.csv", DataFrame(output));
    println("    ---- Finish bus files! ")
end

function write_generator_reactive_output(t,val_opt)
    output = hcat(val_opt.Q_gen_traj)
    filename = string("t=",t,"/gen_reactive_output.csv")
    CSV.write(filename, DataFrame(output));

    # CSV.write("gen_reactive_output.csv", DataFrame(output));
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

function write_RSRV_out(P_rsrv_feedback)
    # println(P_rsrv_feedback)
    CSV.write("P_rsrv_feedback.csv", DataFrame(reshape(P_rsrv_feedback,
        length(P_rsrv_feedback),1),
        [:P_rsrv_feedback]));
    println("    ---- Finish P_rsrv_feedback writting files! ")
end

function write_P_gen_out(P_gen)
    # println(P_rsrv_feedback)
    P_gen_sum = sum(P_gen, dims=1)
    println(size(P_gen_sum))
    CSV.write("P_gen.csv", DataFrame(reshape(P_gen_sum,
        length(P_gen_sum),1),
        [:P_gen_sum]));
    println("    ---- Finish P_gen writting files! ")
end

function write_P_bus_out(P_bus)
    # println(P_rsrv_feedback)
    P_bus_sum = sum(P_bus, dims=1)
    println(size(P_bus_sum))
    CSV.write("P_bus.csv", DataFrame(reshape(P_bus_sum,
        length(P_bus_sum),1),
        [:P_bus_sum]));
    println("    ---- Finish P_bus writting files! ")
end

function write_v_bus_out(t, val_opt)
    # println(P_rsrv_feedback)
    # P_bus_sum = sum(P_bus, dims=1)
    v_feeder = val_opt.v_FOL;
    FOL_pred_length = 12;
    v_feeder_1 = reshape(v_feeder[1,:], FOL_pred_length, 1)
    v_feeder_2 = reshape(v_feeder[2,:], FOL_pred_length, 1)
    v_feeder_3 = reshape(v_feeder[3,:], FOL_pred_length, 1)
    v_FOL = hcat(v_feeder_1, v_feeder_2, v_feeder_3)

    filename = string("./result/v_time_", t, ".csv")
    CSV.write(filename, DataFrame(v_FOL,
        [:feeder_1, :feeder_2, :feeder_3]));
    println("    ---- Finish v_bus writting files! ")
end

function write_output_out_2(val_opt, P_rsrv_feed, val_rp, price, Pd_real, filename)
        # write the solar file
    println("===== GML - Write Output File");
    # name=string("Time", current_time, ".csv");
    cost = sum(val_rp.P_gen)*price.lambda_ct*100;
    time = val_opt.time_solve;
    Pg = sum(val_opt.pg.mu_ct);
    B = sum(val_opt.B);
    R = sum(val_opt.R);
    P_cul = sum(val_opt.P_cul)
    P0 = sum(val_rp.P_gen);
    P_rsrv = P_rsrv_feed;
    status = val_opt.terminate_s;
    Pd = sum(Pd_real)
    cost_agg = val_opt.cost_agg;
    RT_data_feeder=hcat(cost, time, Pg, B, R, Pd, P0, P_rsrv, P_cul, status, cost_agg)
    CSV.write(filename, DataFrame(RT_data_feeder,
        [:cost, :time, :Pg, :B, :R, :Pd_net, :P0, :P_rsrv, :P_cul, :status, :cost_agg]));
    println("    ---- Finish writting files! ")
end

function write_B_out(t, B_feedback)
    B_feedback_temp=reshape(B_feedback, length(B_feedback), 1)
    # println(B_feedback_temp)
    filename = string("time_",t,"_B_feedback.csv")
    CSV.write(filename, DataFrame(B_feedback_temp,
        [:B_feedback]));
    println("    ---- Finish B_feedback writting files! ")
end

function read_RSRV_out()
    data_trace = CSV.File("P_rsrv_feedback.csv") |> DataFrame
    P_rsrv_feedback= collect(data_trace[:,Symbol("P_rsrv_feedback")])
    return reshape(P_rsrv_feedback, 1, length(P_rsrv_feedback))
end

function read_B_out(t)
    filename = string("time_",t,"_B_feedback.csv")
    data_trace = CSV.File(filename) |> DataFrame
    B_feedback = collect(data_trace[:,Symbol("B_feedback")])
    return B_feedback
end

function write_output_out(val_opt, P_rsrv_feed, Pd_bus, filename)
        # write the solar file
    println("===== GML - Write Output File");
    # name=string("Time", current_time, ".csv");
    cost = val_opt.Cost_real;
    time = val_opt.time_solve;
    Pg = sum(val_opt.Pg);
    B = sum(val_opt.B);
    R = sum(val_opt.R);
    P_cul = sum(val_opt.P_cul)
    P0 = sum(val_opt.P_gen);
    P_rsrv = P_rsrv_feed;
    status = val_opt.terminate_s;
    Pd = sum(Pd_bus)
    RT_data_feeder=hcat(cost, time, Pg, B, R, Pd, P0, P_rsrv, P_cul, status)
    CSV.write(filename, DataFrame(RT_data_feeder,
        [:cost, :time, :Pg, :B, :R, :Pd, :P0, :P_rsrv, :P_cul, :status]));
    println("    ---- Finish writting files! ")
end

# function reference_point(Pd_bus, Qd_bus)
#     silence()
#     # network_data = PowerModels.parse_matpower("data/case5.m")
#     network_data = PowerModels.parse_matpower("GML_data/NYISO-data/case_nyiso.m")

#     # Change data

#     # println(network_data["load"])
#     # this is for PQ bus, in which we change the demand but keep the solar;
#     for (load_name, data) in network_data["load"]
#         load_id = parse(Int64, load_name);
#         load_in_bus = data["load_bus"];
#         data["qd"] = Qd_bus[load_in_bus,1];
#         data["pd"] = Pd_bus[load_in_bus,1];
#     end



#     # pm = instantiate_model(network_data, ACPPowerModel, PowerModels.build_opf)
#     # result = optimize_model!(pm, optimizer=Ipopt.Optimizer)
#     result = run_opf(network_data, ACPPowerModel, Ipopt.Optimizer);
#     # admittance = calc_admittance_matrix(network_data)

#     NoBranch = length(result["solution"]["branch"]);
#     NoBus = length(result["solution"]["bus"]);

#     #
#     v_bus = zeros(ComplexF64, NoBus, 1)
#     V_bus = zeros(NoBus, 1)
#     for (bus_name, data) in result["solution"]["bus"]
#         bus_id = parse(Int64, bus_name);
#         # println(bus_id)
#         v_magnitude = result["solution"]["bus"][bus_name]["vm"];
#         # println(v_magnitude)
#         v_polar = result["solution"]["bus"][bus_name]["va"];
#         # println(v_polar)
#         v_bus[bus_id,1] = v_magnitude*exp(v_polar*im);
#         # println(v_bus[bus_id,1])
#         V_bus[bus_id,1] = v_bus[bus_id,1]*conj(v_bus[bus_id,1]);
#         # println(V_bus[bus_id,1])
#     end



#     P_branch=zeros(NoBranch,1);
#     Q_branch=zeros(NoBranch,1);
#     i_branch=zeros(ComplexF64, NoBranch,1);
#     L_branch=zeros(NoBranch,1);

#     for (branch_name, data) in result["solution"]["branch"]
#         # println(branch_name)
#         branch_id = parse(Int64, branch_name);
#         f_bus = network_data["branch"][branch_name]["f_bus"];
#         t_bus = network_data["branch"][branch_name]["t_bus"];
#         P_branch[branch_id,1] = result["solution"]["branch"][branch_name]["pf"];
#         Q_branch[branch_id,1] = result["solution"]["branch"][branch_name]["qf"];
#         i_branch[branch_id,1] = conj(
#         (P_branch[branch_id,1]+Q_branch[branch_id,1]*im)/v_bus[f_bus,1]
#         );
#         L_branch[branch_id,1] = i_branch[branch_id,1]*conj(i_branch[branch_id,1]);

#     end
#     return optimal_setpoints=(p_star=(P_branch), q_star=(Q_branch), l_star=(L_branch),
#         v_star=(V_bus));
# end

function reference_point(Pd_sum, Pd_bus, Qd_bus)
    silence()
    P_Percent = [0.181935169, 0.002855468, 0.001789426, 0.004736269, 0.030496394,
        0.002840238, 0.000639625, 0.386858757, 0.387848653];
    # Pd_sum = sum(Pd_bus);
    # network_data = PowerModels.parse_matpower("data/case5.m")
    network_data = PowerModels.parse_matpower("GML_data/NYISO-data/case_nyiso.m")

    # Change data

    # println(network_data["load"])
    # this is for PQ bus, in which we change the demand but keep the solar;
    for (load_name, data) in network_data["load"]
        load_id = parse(Int64, load_name);
        load_in_bus = data["load_bus"];
        data["qd"] = Qd_bus[load_in_bus,1];
        data["pd"] = Pd_bus[load_in_bus,1];
    end

    for (gen_name, gen_data) in network_data["gen"]
        gen_id = parse(Int64, gen_name);
        gen_data["pg"] = Pd_sum*P_Percent[gen_id];
    end

    # result = run_pf(network_data, ACPPowerModel, Ipopt.Optimizer);
    power_flow_result = run_pf(network_data, ACPPowerModel, Ipopt.Optimizer);
    update_data!(network_data, power_flow_result["solution"])
    flows = calc_branch_flow_ac(network_data)
    update_data!(network_data, flows)

    NoBranch = length(flows["branch"]);
    NoBus = length(power_flow_result["solution"]["bus"]);

    #
    v_bus = zeros(ComplexF64, NoBus, 1)
    V_bus = zeros(NoBus, 1)
    for (bus_name, data) in power_flow_result["solution"]["bus"]
        bus_id = parse(Int64, bus_name);
        # println(bus_id)
        v_magnitude = power_flow_result["solution"]["bus"][bus_name]["vm"];
        # println(v_magnitude)
        v_polar = power_flow_result["solution"]["bus"][bus_name]["va"];
        # println(v_polar)
        v_bus[bus_id,1] = v_magnitude*exp(v_polar*im);
        # println(v_bus[bus_id,1])
        V_bus[bus_id,1] = v_bus[bus_id,1]*conj(v_bus[bus_id,1]);
        # println(V_bus[bus_id,1])
    end

    P_branch=zeros(NoBranch,1);
    Q_branch=zeros(NoBranch,1);
    i_branch=zeros(ComplexF64, NoBranch,1);
    L_branch=zeros(NoBranch,1);

    for (branch_name, data) in network_data["branch"]
        # println(branch_name)
        branch_id = parse(Int64, branch_name);
        f_bus = network_data["branch"][branch_name]["f_bus"];
        t_bus = network_data["branch"][branch_name]["t_bus"];
        P_branch[branch_id,1] = network_data["branch"][branch_name]["pf"];
        Q_branch[branch_id,1] = network_data["branch"][branch_name]["qf"];
        i_branch[branch_id,1] = conj(
        (P_branch[branch_id,1]+Q_branch[branch_id,1]*im)/v_bus[f_bus,1]
        );
        L_branch[branch_id,1] = i_branch[branch_id,1]*conj(i_branch[branch_id,1]);

    end

    P_gen=zeros(9,1);
    for (gen_name, data) in power_flow_result["solution"]["gen"]
        gen_id = parse(Int64, gen_name);
        P_gen[gen_id, 1] = power_flow_result["solution"]["gen"][gen_name]["pg"]
    end

    return optimal_setpoints=(p_star=(P_branch), q_star=(Q_branch), l_star=(L_branch),
        v_star=(V_bus), power_flow_result=(power_flow_result),
        network_data=(network_data), P_gen=(P_gen));
end
