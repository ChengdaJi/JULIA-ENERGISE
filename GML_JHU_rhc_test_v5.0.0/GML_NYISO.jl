function GML_Sys_Ava_NYISO(T, pd, ancillary_type, B_cap, icdf, bus_struct,
    shunt_struct, branch_struct, gen_struct, baseMVA)
    println("===== GML - Boundaries Buildup");
    ###############################################################################
    # shunt level
    T=288;
    NoBus = length(bus_struct.baseKV)
    Qf_max = zeros(NoBus,T)
    Qf_min = zeros(NoBus,T)
    for bus = 1:NoBus
        # println(bus)
        if bus_struct.type[bus] == 1
            Qf_max[bus, :]=1.5*positive_array(pd.traj[bus, :]);
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
    # B_max = zeros(NoShunt, T)
    B_min = zeros(NoBus, T);
    R_rate = 1/3;
    R_max = R_rate*B_max/baseMVA;
    R_min = -R_max;
    # println(R_max[:,1])
    # println(R_min[:,1])
    W = zeros(NoShunt, T);
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
    println("===== GML - finish building boundaries =====")
    obj=(T=(T), icdf = (icdf), Qf_max=(Qf_max), Qf_min=(Qf_min), Pg_min=(Pg_min),
    B_max=(B_max), B_min=(B_min), R_max=(R_max), R_min=(R_min),
    tau=(tau), P_rsrv_min=(P_rsrv_min), k=(k),V_max=(V_max), V_min=(V_min),
    r=(r), x=(x), beta=(beta));

    return obj
end

function optimal_NYISO(NoShunt, SN, current_time, obj, ancillary_type,
    feedback, pd, pg, price, shunt_data, bus_data, branch_data, gen_data);

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

    NoBus=length(bus_data.baseKV)
    BrN=length(branch_data.fbus)
    # NoGen=length(gen_data.id)
    NoGen=9;

    m = Model(with_optimizer(Mosek.Optimizer, QUIET=true, MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-3,
    MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-3))
    # define the real-time variables

    # Shunt level
    @variable(m, Pg_rt[1:NoShunt, 1])
    @variable(m, Qf_rt[1:NoShunt, 1])
    @variable(m, B_rt[1:NoShunt, 1])
    @variable(m, R_rt[1:NoShunt, 1])
    # Bus level
    @variable(m, P_bus_rt[1:NoBus, 1])
    @variable(m, Q_bus_rt[1:NoBus, 1])
    @variable(m, v_rt[1:NoBus, 1])
    # Branch level
    @variable(m, P_br_rt[1:BrN, 1])
    @variable(m, Q_br_rt[1:BrN, 1])

    @variable(m, P_gen_rt[1:NoGen, 1])
    @variable(m, Q_gen_rt[1:NoGen, 1])


    # println(" ---- Real Time Constraint Buildup ")
    for shunt=1:NoShunt
        if shunt_data.type[shunt] == 1
            shuntbMVA=shunt_struct.baseMVA[shunt]
            @constraint(m, 0<=Pg_rt[shunt,1]);
            @constraint(m, Pg_rt[shunt,1]<=
                positive_scalar(icdf*sqrt(pd.sigma[shunt,1]+pg.sigma[shunt,1])/shuntbMVA
                +pg.mu[shunt,1]/shuntbMVA));
            @constraint(m, R_min[shunt,1]<= R_rt[shunt,1]);
            @constraint(m, R_rt[shunt,1]<= R_max[shunt,1]);
            @constraint(m, B_rt[shunt,1]==B_feedback[shunt,1]);
        else
            @constraint(m, Pg_rt[shunt,1]==0);
            @constraint(m, R_rt[shunt,1]==0);
            @constraint(m, B_rt[shunt,1]==0);
        end
        @constraint(m, Qf_min[shunt,1]<=Qf_rt[shunt,1]);
        @constraint(m, Qf_rt[shunt,1]<=Qf_max[shunt,1]);
    end


    for Bus=1:NoBus
        @constraint(m, V_min^2<= v_rt[Bus,1]);
        @constraint(m, v_rt[Bus,1]<= V_max^2);
        if bus_data.type[Bus]==1
            bus_shunt_list =
            setdiff(findall(hasshunt->hasshunt==Bus, shunt.find_bus),198:209);
            if ~isempty(bus_shunt_list)
                @constraint(m, P_bus_rt[Bus,1]==sum(
                Pd[shunt,1]/shunt_struct.baseMVA[shunt]
                -Pg_rt[shunt,1]-R_rt[shunt,1] for shunt in bus_shunt_list));
            else
                @constraint(m, P_bus_rt[Bus,1]==0)
            end
        else
            @constraint(m, P_bus_rt[Bus,1]==0)
        end
    end

    for Brh=1:BrN
        fbus = branch_data.fbus[Brh];
        tbus = branch_data.tbus[Brh];
        tbus_out_branch_list =
            setdiff(findall(isafbus->isafbus==tbus, branch_data.fbus), Brh);
        tbus_in_branch_list =
            setdiff(findall(isafbus->isafbus==tbus, branch_data.tbus), Brh);
        if tbus in gen_struct.bus
            genid_list = findall(isagen->isagen==tbus, gen_data.bus)
            genid=genid_list[1]
            if ~isempty(tbus_out_branch_list) && ~isempty(tbus_in_branch_list)
                @constraint(m,P_br_rt[Brh,1]==
                    -P_gen_rt[genid, 1]
                    +sum(P_br_rt[outbrunch, 1] for outbrunch in tbus_out_branch_list)
                    -sum(P_br_rt[inbrunch, 1] for inbrunch in tbus_in_branch_list));

            elseif isempty(tbus_out_branch_list) && ~isempty(tbus_in_branch_list)
                @constraint(m,P_br_rt[Brh,1]==
                    -P_gen_rt[genid, 1]
                    -sum(P_br_rt[inbrunch, 1] for inbrunch in tbus_in_branch_list));

            elseif ~isempty(tbus_out_branch_list) && isempty(tbus_in_branch_list)
                @constraint(m,P_br_rt[Brh,1]==
                    -P_gen_rt[genid, 1]
                    +sum(P_br_rt[outbrunch, 1] for outbrunch in tbus_out_branch_list));

            elseif isempty(tbus_out_branch_list) && isempty(tbus_in_branch_list)
                @constraint(m,P_br_rt[Brh,1]==
                    -P_gen_rt[genid, 1])
            end
        else
            if ~isempty(tbus_out_branch_list) && ~isempty(tbus_in_branch_list)
                @constraint(m,P_br_rt[Brh,1]==
                    P_bus_rt[tbus, 1]
                    +sum(P_br_rt[outbrunch,1] for outbrunch in tbus_out_branch_list)
                    -sum(P_br_rt[inbrunch,1] for inbrunch in tbus_in_branch_list));

            elseif isempty(tbus_out_branch_list) && ~isempty(tbus_in_branch_list)
                @constraint(m,P_br_rt[Brh,1]==
                    P_bus_rt[tbus, 1]
                    -sum(P_br_rt[inbrunch,1] for inbrunch in tbus_in_branch_list));

            elseif ~isempty(tbus_out_branch_list) && isempty(tbus_in_branch_list)
                @constraint(m,P_br_rt[Brh,1]==
                    P_bus_rt[tbus, 1]
                    +sum(P_br_rt[outbrunch,1] for outbrunch in tbus_out_branch_list));

            elseif isempty(tbus_out_branch_list) && isempty(tbus_in_branch_list)
                @constraint(m,P_br_rt[Brh,1]==
                    P_bus_rt[tbus, 1]);
            end
        end

        # @constraint(m,v_rt[tbus,1]==
        #     v_rt[fbus,1]-2*(r[Brh]*P_br_rt[Brh,1]+x[Brh]*Q_br_rt[Brh,1]));

    end


    for Gen=1:NoGen
        gen_out_branch_list=findall(isafbus->isafbus==Gen, branch_data.fbus);
        gen_in_branch_list=findall(isatbus->isatbus==Gen, branch_data.tbus);
        if ~isempty(gen_out_branch_list) && ~isempty(gen_in_branch_list)
        @constraint(m, P_gen_rt[Gen,1]==
           sum(P_br_rt[outbranch,1] for outbranch in gen_out_branch_list)
           -sum(P_br_rt[inbranch,1] for inbranch in gen_in_branch_list));
        elseif isempty(gen_out_branch_list) && ~isempty(gen_in_branch_list)
            @constraint(m, P_gen_rt[Gen,1]==
            -sum(P_br_rt[inbranch,1] for inbranch in gen_in_branch_list));
        elseif ~isempty(gen_out_branch_list) && isempty(gen_in_branch_list)
            @constraint(m, P_gen_rt[Gen,1]==
            sum(P_br_rt[outbranch,1] for outbranch in gen_out_branch_list));
        # elseif isempty(gen_out_branch_list) && isempty(gen_in_branch_list)
        #     @constraint(m, P_gen_rt[Gen,1]== 0);
        #     println("???")
        end

        @constraint(m,
            P_gen_rt[Gen,1]<=gen_data.Pmax[Gen]/bus_data.baseMVA[Gen])
        @constraint(m,
            P_gen_rt[Gen,1]>=-gen_data.Pmax[Gen]/bus_data.baseMVA[Gen])
        # @constraint(m,
        #     Q_gen_rt[Gen,1]<=gen_data.Qmax[Gen]/bus_data.baseMVA[Gen])
        # @constraint(m,
        #     Q_gen_rt[Gen,1]>=gen_data.Qmin[Gen]/bus_data.baseMVA[Gen])
    end
    # shunt
    @variable(m, Pg[1:SN, 1:NoShunt, 1:T-1]); # the real power output
    @variable(m, Qf[1:SN, 1:NoShunt, 1:T-1]); # the real power output
    @variable(m, B[1:SN, 1:NoShunt, 1:T-1]); # the storage
    @variable(m, R[1:SN, 1:NoShunt, 1:T-1]);# the charge/discharge rate
    # bus
    @variable(m, P_bus[1:SN, 1:NoBus, 1:T-1])
    @variable(m, Q_bus[1:SN, 1:NoBus, 1:T-1])
    @variable(m, v[1:SN, 1:NoBus, 1:T-1])
    # Branch level
    @variable(m, P_br[1:SN,1:BrN, 1:T-1])
    @variable(m, Q_br[1:SN,1:BrN, 1:T-1])

    @variable(m, P_gen[1:SN, 1:NoGen, 1:T-1])
    @variable(m, Q_gen[1:SN, 1:NoGen, 1:T-1])
    # #

    if ancillary_type == "10min" || ancillary_type == "30min"
        @objective(m, Min,
            fn_cost_RHC_anc(delta_t,P_bus_rt,P_bus,Pg_rt,Pg,P_rsrv_rt,
            P_rsrv,price,pg,pd,beta,SN,obj, gen_struct, bus_struct, branch_struct))
    else
        @objective(m, Min,
            fn_cost_RHC_rt(delta_t, P_gen, P_gen_rt, Pg_rt,Pg,price,
                pg,pd,beta,SN,obj, gen_struct, bus_struct, branch_struct))
    end
    #
    status=optimize!(m);
    #
    println(string("    ----", termination_status(m)))
    # println(MOI.PrimalStatus())
    # println(MOI.DualStatus())
    cost_o = JuMP.objective_value(m);
    time_solve=MOI.get(m, MOI.SolveTime());
    println(string("    ----Solve Time: ", time_solve))
    println(string("    ----Optimal Cost: ", cost_o))
    # ## obtaining value
    # Pg_rt_o=JuMP.value.(Pg_rt)
    # Pg_rt_s=JuMP.value.(Pg)
    # #
    # # R_rt_o=JuMP.value.(R_rt)
    # # R_s=JuMP.value.(R)
    # # R_total = hcat(R_rt_o, R_s[1,:,:])
    # #
    # P_bus_ct_o=JuMP.value.(P_bus_rt)
    # P_bus_sum_ct_o = sum(P_bus_ct_o[i,1] for i in 1:79);
    # P_bus_o=JuMP.value.(P_bus)
    # gen_list = gen_struct.id;
    # for gen=1:length(gen_list)
    #     println(gen_list[gen])
    #     println(P_bus_rt_o[gen_list[gen], 1])
    #     println(P_bus_rt_o[gen_list[gen], 1]*baseMVA[gen_list[gen]])
    #     println(baseMVA[gen_list[gen]])
    # end
    P_gen_rt_o=JuMP.value.(P_gen_rt)
    println(sum(P_gen_rt_o))
    println(P_gen_rt_o)
    P_br_rt_o=JuMP.value.(P_br_rt)
    println(sum(P_br_rt_o))


    for Gen=6:6
        gen_out_branch_list=findall(isafbus->isafbus==Gen, branch_data.fbus);
        gen_in_branch_list=findall(isafbus->isafbus==Gen, branch_data.tbus);
        println(string("out branch", gen_out_branch_list))
        println(string("in branch", gen_in_branch_list))
        println(sum(P_br_rt_o[outbranch,1] for outbranch in gen_out_branch_list));
        # println(sum(P_br_rt_o[inbranch,1] for inbranch in gen_in_branch_list));
    end
    shuntbMVA = shunt_struct.baseMVA;
    Pd_pu=zeros(209,1)
    for Shunt = 1:197
        Pd_pu[Shunt,1]=Pd[Shunt,1]/shuntbMVA[Shunt]
    end
    println(sum(Pd_pu))

    Pg_rt_o=JuMP.value.(Pg_rt);
    println(sum(Pg_rt_o))
    Pg_bus_rt_o=JuMP.value.(P_bus_rt);
    println(sum(Pg_bus_rt_o))
    # P_gen_o=JuMP.value.(P_gen)
    # P_gen_traj=hcat(P_gen_rt_o, P_gen_o[1,:,:]);
    # println(sum(P_gen_traj, dims=1))
    # println(P_gen_rt_o)
    # println(size(P_br_rt))
    # P_br_rt_o=JuMP.value.(P_br_rt)
    # println(size(P_br_rt_o))
    # println(sum(P_br_rt_o))
    println(sum(Pd[1:197,1]))
    # P_gen_sum_o = reshape(
    #     sum(P_bus_o[1,gen, :]*baseMVA[gen] for gen in gen_list),T-1,1)
    # P0 = vcat(P_gen_sum_ct_o, P_gen_sum_o)
    # Pg_rt_o=JuMP.value.(Pg_rt)
    # Pg_o=JuMP.value.(Pg)
    # Pg_sum_pu = hcat(Pg_rt_o,Pg_o[1,:,:])
    # Pg_traj = zeros(209,1)
    # for shunt =1:209
    #     Pg_traj[shunt,1] = Pg_rt_o[shunt,1]*shunt_struct.baseMVA[shunt]
    # end
    # println(Pg_sum)
    # println(sum(Pg_traj)/sum(Pd[:,1]))
    # println(sum(P_gen_sum_ct_o)/sum(Pd[:,1]))
    # println(sum(P_gen_sum_ct_o_2)/sum(Pd[:,1]))
    #
    # # P_bus_3 = hcat(P_bus_rt_o, P_bus_s[3,:,:])
    # # println(sum(P_bus_3[1,:]))
    # # P_bus_4 = hcat(P_bus_rt_o, P_bus_s[4,:,:])
    # # println(sum(P_bus_4[1,:]))
    # # P_bus_5 = hcat(P_bus_rt_o, P_bus_s[5,:,:])
    # # println(sum(P_bus_5[1,:]))
    # # P_bus_6 = hcat(P_bus_rt_o, P_bus_s[6,:,:])
    # # println(sum(P_bus_6[1,:]))
    #
    # # P_bus_s=JuMP.value.(P_bus)
    # # P_bus_total = hcat(P_bus_rt_o, P_bus_s[1,:,:])
    # #
    # # if ancillary_type == "10min" || ancillary_type == "30min"
    # #     P_rsrv_o=JuMP.value(P_rsrv_rt);
    # #     P_rsrv_s=JuMP.value.(P_rsrv);
    # #     P_rsrv_total = hcat([P_rsrv_o], P_rsrv_s);
    # # else
    # #     P_rsrv_total = zeros(1,T);
    # # end
    #
    # # val_opt = (lambda1=(price.lambda_scenario[1, 1:T]/12),P_gen1=(P_bus_1[1,:]*100),
    # # lambda2=(price.lambda_scenario[2, 1:T]/12), P_gen2=(P_bus_2[1,:]*100))
    # # val_opt = (P_bus = (P_bus_total), Pg = (Pg_total), R=(R_total), P_rsrv = (P_rsrv_total),
    # #     Cost = (cost_o), Solve_time=(time_solve))


    return 0
end


function fn_cost_RHC_rt(delta_t, P_gen, P_gen_rt, Pg_rt,Pg,price,
    pg,pd,beta,SN,obj, gen_struct, bus_struct, branch_struct)

    T=obj.T;
    icdf = obj.icdf;

    baseMVA=gen_struct.baseMVA;
    sum_prob = sum(price.probability[1:SN])
    gen_list = gen_struct.id;

    P_gen_sum_ct = sum(P_gen_rt[gen, 1]*baseMVA[gen] for gen in gen_list)
    # P_gen_sum_ct = sum(P_gen_rt[:, 1])
    lambda_ct = price.probability[1]/sum_prob*price.lambda_scenario[1,1]
    #
    # Pg_ct_sum_diff = sum(Pg_rt) -
    #     sum(positive_array(icdf.*sqrt.(pd.sigma[:,1]+pg.sigma[:,1]))
    #     +pg.mu[:,1])
    # Cost_Pg_ct_diff = Pg_ct_sum_diff*beta*delta_t;
    #
    # if SN>1
    #     for scenario=2:SN
    #         lambda_ct =lambda_ct+
    #         price.probability[scenario]/sum_prob*price.lambda_scenario[scenario,1];
    #     end
    # end
    Cost_P_gen_sum_ct = delta_t*lambda_ct*P_gen_sum_ct;
    # #
    # # # =======
    P_gen_sum_scenario = reshape(
        sum(P_gen[1, gen, :]*baseMVA[gen] for gen in gen_list),T-1,1)
    # # #
    lambda_scenario = reshape(
        price.probability[1]/sum_prob*price.lambda_scenario[1,2:T],1,T-1)
    # # #
    Cost_P_gen_sum_scenario = delta_t* lambda_scenario * P_gen_sum_scenario;
    # #
    # # Pg_prob_sum_scenario = price.probability[1]/sum_prob*sum(Pg[1, :, :])
    # #
    # # if SN>1
    # #     for scenario=2:SN
    # #         P_gen_sum_scenario = reshape(
    # #         sum(P_bus[scenario, gen, :]*baseMVA[gen] for gen in gen_list),T-1,1)
    # #
    # #         lambda_scenario_prob = reshape(
    # #         price.probability[scenario]/sum_prob*price.lambda_scenario[scenario,2:T],1,T-1)
    # #
    # #         Cost_P_gen_sum_scenario = Cost_P_gen_sum_scenario+delta_t* lambda_scenario_prob * P_gen_sum_scenario;
    # #
    # #         # Pg_prob_sum_scenario = Pg_prob_sum_scenario+price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :])
    # #     end
    # # end
    # # Pg_diff_sum_scenario = Pg_prob_sum_scenario
    # #     -positive_array(icdf.*sqrt.(pd.sigma[:,:]+pg.sigma[:,:])+pg.mu[:,:])
    # #
    # # Cost_Pg_scenario_diff = Pg_diff_sum_scenario*beta*delta_t;
    # #
    # #
    # # println(string("    ----Case: Real-time Balancing"))
    # # return (Cost_P_gen_sum_ct+Cost_P_gen_sum_scenario[1,1]
    # # +Cost_Pg_ct_diff+Cost_Pg_scenario_diff)*base.MVA
    return Cost_P_gen_sum_ct
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

function read_NY_demand_data(shunt, bus)
    P=shunt.P;
    NoBus = length(bus.baseKV)


    demand_unit=matread("../data/NYISO-data/normalized_demand.mat")["normalized_demand"]
    demand=zeros(length(P),576)
    for n_shunt=1:length(P)
        demand[n_shunt,:] = reshape(
        demand_unit[:, n_shunt].*(576).*P[n_shunt].*0.7, 1,576);
    end

    demand_da_unit=matread("../data/NYISO-data/normalized_demand_da.mat")["normalized_demand_da"]
    demand_da=zeros(length(P),576)
    for n_shunt=1:length(P)
        demand_da[n_shunt,:] = reshape(
        demand_da_unit[:, n_shunt].*(576).*P[n_shunt].*0.7, 1,576);
    end

    demand_bus = zeros(NoBus, 576)
    demand_da_bus = zeros(NoBus, 576)

    for bus = 1:NoBus
        shunt_list = findall(inbus->inbus==bus, shunt.find_bus)
        if isempty(shunt_list)
            demand_bus[bus, :]=zeros(1, 576)
            demand_da_bus[bus, :]=zeros(1, 576)
        else
            demand_bus[bus, :] = sum(demand[shunt, :] for shunt in shunt_list)
            demand_da_bus[bus, :] = sum(demand_da[shunt, :] for shunt in shunt_list)
        end
    end

    pd_raw = (pd_rt = (demand_bus), pd_da = (demand_da_bus));
    return pd_raw
end

function read_NY_solar_data(shunt, bus)
    P=shunt.P;
    NoBus = length(bus.baseKV);


    solar_unit=matread("../data/NYISO-data/normalized_solar.mat")["normalized_solar"]
    solar=zeros(length(P),576)
    for n_shunt=1:length(P)
        solar[n_shunt,:] = reshape(
        solar_unit[:, n_shunt].*(576).*P[n_shunt].*0.7, 1,576);
    end

    solar_da_unit=matread("../data/NYISO-data/normalized_solar_da.mat")["normalized_solar_da"]
    solar_da=zeros(length(P),576)
    for n_shunt=1:length(P)
        solar_da[n_shunt,:] = reshape(
        solar_da_unit[:, n_shunt].*(576).*P[n_shunt].*0.7, 1,576);
    end
    solar_bus = zeros(NoBus, 576)
    solar_da_bus = zeros(NoBus, 576)
    for bus = 1:NoBus
        shunt_list = findall(inbus->inbus==bus, shunt.find_bus)
        if isempty(shunt_list)
            solar_bus[bus, :]=zeros(1, 576)
            solar_da_bus[bus, :]=zeros(1, 576)
        else
            solar_bus[bus, :] = sum(solar[shunt, :] for shunt in shunt_list)
            solar_da_bus[bus, :] = sum(solar_da[shunt, :] for shunt in shunt_list)
        end
    end
    pg_raw = (pg_rt = (solar_bus), pg_da = (solar_da_bus));
    return pg_raw
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
