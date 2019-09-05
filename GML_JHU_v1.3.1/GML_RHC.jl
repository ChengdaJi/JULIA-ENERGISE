function GML_RHC_DET(current_time, B)
    B_feedback = repeat([B], 12,1)
    ancillary_type="without"
    # # of timeslots
    T=288;
    # BN=252;
    BN=3;
    # # of feeders
    F=BN*4;
    # # of secnarios
    SN=1;
    # penetation level
    p_rate=0.5; # 0.25 0.5 0.75 1
    # # prediction Pred_length
    Pred_length=24; # 5 min per slots
    # load price data
    price_raw = read_price_data()
    # delta_rt_raw=matread("../data/price_prediction.mat");
    ## read read_demand_data
    pd_raw = read_demand_data()
    # pd_noise = matread("../data/demand_noise.mat")["demand_noise"];
    pg_raw = read_solar_data()
    # # initial SOC
    P_rsrv_feedback = [];
    # #
    ct_printout = string("    GML - At Time ", current_time);
    println("=================================================")
    println(ct_printout)
    feedback = feedback_struct(B_feedback, [])
    price = price_traj_det(current_time, ancillary_type, price_raw, T)
    pd = pd_traj(current_time, pd_raw, T);
    pg = pg_traj(current_time, pg_raw, T);
    obj = GML_Sys_Ava(T, F, BN, SN, pd, ancillary_type, 0);
    val_opt = optimal_stoach_scenario_offline(current_time, obj, feedback, pd, pg,
    price, ancillary_type);
    println("=================================================")
    return val_opt
end


function GML_Sys_Ava(T, F, BN, SN, pd, ancillary_type, icdf)

    # println("    GML - Boundaries Buildup");
    ###############################################################################

    # feeder level
    Qf_max=0.05*positive_array(pd.traj);
    Qf_min = -Qf_max;
    # minimum solar
    Pg_min = zeros(F, T);
    # battery
    B_max = ones(F,1)*ones(1,T)*0.5;
    B_min = zeros(F,T);
    R_max = 2.325*ones(F,T);
    R_min = -R_max;
    W=zeros(F,T);
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
    S=35*ones(1,T);
    # V
    Base_V=69;
    V_min = Base_V*0.95;
    V_max = Base_V*1.05;

    # impedance
    r=[0.13275; 0.13275; 0.199125];
    x=[1.00426; 1.00426; 1.50639];
    # penalty
    beta=-15;
    # println("===== GML - finish building boundaries =====")
    obj=obj_struct(T,F,BN,SN, Pg_min,Qf_max,Qf_min,B_max,B_min,R_max,R_min,W,
    delta_t,P_rsrv_min,tau,k,V_max,V_min,r,x,beta,S,icdf);
    return obj
end

function optimal_stoach_scenario_offline(current_time, obj, feedback, pd, pg,
price, ancillary_type);
    # println("    GML - Optimization")

    #################
    #  F,T or BN, T #
    #################
    ## The parameters

    T = obj.T; # number of time slots
    F = obj.F; # number of feeders
    BN = obj.BN; # number of banks
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

    TB=[];
    for tb_index=1:BN
        push!(TB, (tb_index-1)*4+1:tb_index*4)
    end
    V_min_0=69*0.95;
    V_max_0=69*1.05;



    m = Model(with_optimizer(Mosek.Optimizer, QUIET=true, MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-3,
    MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-3))
    # define the real-time variables
    # Feeder level
    @variable(m, Pg_rt[1:F, 1])
    @variable(m, Qf_rt[1:F, 1])
    @variable(m, B_rt[1:F, 1])
    @variable(m, R_rt[1:F, 1])
    # Bank level
    @variable(m, P_hat_rt[1:BN, 1])
    @variable(m, Q_hat_rt[1:BN, 1])
    @variable(m, l_rt[1:BN, 1])
    @variable(m, v_rt[1:BN, 1])
    @variable(m, v_0_rt[1,1])
    @variable(m, P_rt[1:BN, 1])
    @variable(m, Q_rt[1:BN, 1])

    # @variable(m, u_temp_rt[1:BN, 1])
    # @variable(m, v_temp_rt[1:BN, 1])

    # println(" ---- Real Time Constraint Buildup ")
    for feeder=1:F
        @constraint(m, 0<=Pg_rt[feeder,1]);
        @constraint(m, pg.mu_rt[feeder,1]>=Pg_rt[feeder,1]);
        @constraint(m, R_min[feeder,1]<= R_rt[feeder,1]);
        @constraint(m, R_rt[feeder,1]<= R_max[feeder,1]);
        @constraint(m, B_rt[feeder,1]==B_feedback[feeder,1]);
        @constraint(m, Qf_min[feeder,1]<=Qf_rt[feeder,1]);
        @constraint(m, Qf_rt[feeder,1]<=Qf_max[feeder,1]);
    end

    for bank=1:BN
        @constraint(m, P_rt[bank,1]==
            sum(Pd[TB[bank][1]:TB[bank][end],1])-
            sum(Pg_rt[TB[bank][feeder],1]-R_rt[TB[bank][feeder],1] for feeder=1:4));
        @constraint(m, Q_rt[bank,1]==sum(Qf_rt[TB[bank][feeder_ite],1] for feeder_ite=1:4));
        @constraint(m, [0.5*S[1,1], S[1,1], P_rt[bank,1], Q_rt[bank,1]] in RotatedSecondOrderCone());
        @constraint(m,P_hat_rt[bank,1]==
            r[bank]*l_rt[bank,1]+P_rt[bank,1]);
        @constraint(m,Q_hat_rt[bank,1]==
            x[bank]*l_rt[bank,1]+Q_rt[bank,1]);
        @constraint(m,v_rt[bank,1]==v_0_rt[1,1]+
            (r[bank]^2+x[bank]^2)*l_rt[bank,1]-
            2*(r[bank]*P_hat_rt[bank,1]+x[bank]*Q_hat_rt[bank,1]));
        @constraint(m, [0.5*l_rt[bank,1], v_0_rt[1,1], P_hat_rt[bank,1], Q_hat_rt[bank,1]] in RotatedSecondOrderCone());
        @constraint(m, V_min^2<= v_rt[bank,1]);
        @constraint(m, v_rt[bank,1]<= V_max^2);
    end
    @constraint(m, v_0_rt[1,1]>=V_min^2);
    @constraint(m, v_0_rt[1,1]<=V_max^2);
    @constraint(m, sum(Q_hat_rt)<=0.05*sum(P_hat_rt));



    # println(" ---- Stochastic Constraint Buildup")
    @variable(m, Pg[1:F, 1:T-1]); # the real power output
    @variable(m, Qf[1:F, 1:T-1]); # the real power output
    @variable(m, B[1:F, 1:T-1]); # the storage
    @variable(m, R[1:F, 1:T-1]);# the charge/discharge rate
    # bank level
    @variable(m, P_hat[1:BN, 1:T-1]);# the total delivered power on line
    @variable(m, Q_hat[1:BN, 1:T-1]);# the total delivered reactive on line
    @variable(m, l[1:BN, 1:T-1]);# current square
    @variable(m, v[1:BN, 1:T-1]);# voltage
    @variable(m, v_0[1:T-1]);# voltage for P0
    @variable(m, P_[1:BN, 1:T-1]); # power distributed from bank
    @variable(m, Q_[1:BN, 1:T-1]); # reactive power distributed from bank

    for scenario = 1:SN
        for t=1:T-1
            for feeder=1:F
                @constraint(m, Pg_min[feeder, t]<=Pg[feeder ,t]);
                @constraint(m, Pg[feeder ,t]<=pg.mu_scenario[feeder,t]);
                @constraint(m, B_min[feeder,t+1] <= B[feeder,t]);
                @constraint(m, B[feeder,t]<= B_max[feeder,t+1]);
                @constraint(m, R_min[feeder,t+1] <= R[feeder,t]);
                @constraint(m, R[feeder,t]<= R_max[feeder,t+1]);
                if t==1
                    @constraint(m, B[feeder,1] == B_rt[feeder,1]-delta_t*R_rt[feeder,1])
                else
                    @constraint(m, B[feeder,t] ==
                        B[feeder,t-1] - R[feeder,t-1]*delta_t)
                end
                @constraint(m, Qf_min[feeder,t+1] <= Qf[feeder,t]);
                @constraint(m, Qf[feeder,t]<= Qf_max[feeder,t+1]);
            end
            for bank=1:BN
                temp_Pg = Pg[TB[bank][1],t];
                temp_R = R[TB[bank][1],t];
                temp_Qf = Qf[TB[bank][1],t];
                for feeder_ite = 2:length(TB[bank])
                    temp_Pg = temp_Pg +Pg[TB[bank][feeder_ite],t];
                    temp_R = temp_R+R[TB[bank][feeder_ite],t];
                    temp_Qf = temp_Qf+Qf[TB[bank][feeder_ite],t];
                end
                @constraint(m, P_[bank,t]==
                    sum(Pd[TB[bank][1]:TB[bank][end],t])
                    -temp_Pg-temp_R);
                @constraint(m, Q_[bank,t]==temp_Qf);
                @constraint(m, [0.5*S[1,t+1], S[1,t+1], P_[bank,t], Q_[bank,t]] in RotatedSecondOrderCone());
                @constraint(m,P_hat[bank,t]==
                    r[bank]*l[bank,t]+P_[bank,t]);
                @constraint(m,Q_hat[bank,t]==
                    x[bank]*l[bank,t]+Q_[bank,t]);
                @constraint(m,v[bank,t]==v_0[t]+
                    (r[bank]^2+x[bank]^2)*l[bank,t]-
                    2*(r[bank]*P_hat[bank,t]+x[bank]*Q_hat[bank,t]));
                @constraint(m, [0.5*l[bank,t], v_0[t], P_hat[bank,t], Q_hat[bank,t]] in RotatedSecondOrderCone());
                @constraint(m, V_min^2<=v[bank,t]);
                @constraint(m, V_max^2>=v[bank,t]);
            end
            @constraint(m, V_min_0^2<=v_0[t]);
            @constraint(m, V_max_0^2>=v_0[t]);
            @constraint(m, sum(Q_hat[:,t])<=0.05*sum(P_hat[:,t]));
        end
    end
    @objective(m, Min,
        fn_cost_RHC_det(delta_t,price,P_hat_rt,P_hat,Pg_rt,Pg,pg,pd,beta,obj))
    status=optimize!(m);

    # println(string("    ----", termination_status(m)))
    # println(MOI.PrimalStatus())
    # println(MOI.DualStatus())
    cost_o = JuMP.objective_value(m);
    time_solve=MOI.get(m, MOI.SolveTime());
    # println(string("    ----Solve Time: ", time_solve))
    # println(string("    ----Optimal Cost: ", cost_o))
    ## obtaining value
    # real_time
    Qf_o=JuMP.value.(Qf_rt)
    Pg_o=JuMP.value.(Pg_rt)
    B_o=JuMP.value.(B_rt)
    R_o=JuMP.value.(R_rt)
    P_hat_o=JuMP.value.(P_hat_rt)
    Q_hat_o=JuMP.value.(Q_hat_rt)
    l_o=JuMP.value.(l_rt)
    v_o=JuMP.value.(v_rt)

    Qf_o_s=JuMP.value.(Qf)
    Pg_o_s=JuMP.value.(Pg)
    R_o_s=JuMP.value.(R)
    v_o_s=JuMP.value.(v)


    # Pg_s = JuMP.value.(Pg);
    # Qf_s = JuMP.value.(Qf);
    # B_s = JuMP.value.(B);
    # R_s = JuMP.value.(R);
    # P_hat_s = JuMP.value.(P_hat);
    # Q_hat_s = JuMP.value.(Q_hat);
    # l_s = JuMP.value.(l);
    # v_s = JuMP.value.(v);
    # v_0_s = JuMP.value.(v_0);
    # P_s= JuMP.value.(P_);
    # Q_s= JuMP.value.(Q_);


    P_rsrv_o=0;
    B_rsrv_o=0;
    P_0_o=sum(P_hat_o)
    cost_o=JuMP.objective_value(m);
    Pf_o=zeros(F,1)
    for feeder = 1:F
        Pf_o[feeder,1]=Pd[feeder,1]-Pg_o[feeder,1]-R_o[feeder,1]
    end
    FOL_Horizon = 6;
    R_opt = zeros(1,6);
    Pg_opt = zeros(1,6);
    Qf_opt = zeros(1,6);
    v_opt = zeros(1,6);
    println(size(R_o))
    for FOL_num = 1:FOL_Horizon
        if FOL_num ==1
            R_opt[1, FOL_num]=R_o[1,FOL_num];
            Pg_opt[1, FOL_num]=Pg_o[1,FOL_num];
            Qf_opt[1, FOL_num]=Qf_o[1,FOL_num];
            v_opt[1, FOL_num]=sqrt(v_o[1,FOL_num]);
        else
            R_opt[1, FOL_num]=R_o_s[1,FOL_num];
            Pg_opt[1, FOL_num]=Pg_o_s[1,FOL_num];
            Qf_opt[1, FOL_num]=Qf_o_s[1,FOL_num];
            v_opt[1, FOL_num]=sqrt(v_o_s[1,FOL_num]);
        end
    end
    val_opt = (R=(R_opt), Pg=(Pg_opt), Qf=(Qf_opt), v=(v_opt))
    # val_opt = Optimization_output_struct(Pf_o, Qf_o, Pg_o, B_o, R_o, P_hat_o, Q_hat_o, l_o, v_o, P_rsrv_o, B_rsrv_o, P_0_o, cost_o, time_solve)
    return val_opt
end

function fn_cost_RHC_det(delta_t,price,P_hat_rt,P_hat,Pg_rt,Pg,pg,pd,beta,obj)
    # println("    ---- Load the optimal function")
    T=obj.T;
    P_hat_scenario=sum(P_hat[:, :], dims=1);
    Cost_P_hat_scenario = P_hat_scenario*reshape(price.lambda_scenario[1, :],T-1,1);
    Pg_diff_scenario=sum(Pg_rt[:,1])-sum(pg.mu_rt)+sum(Pg[:, :]);
    Pg_diff_scenario = Pg_diff_scenario-
        sum(pg.mu_scenario);
    return delta_t*price.lambda_rt*sum(P_hat_rt[:,1])+(delta_t*Cost_P_hat_scenario)[1,1]+
        delta_t*beta*Pg_diff_scenario
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

function write_output_out(val_opt, current_time)
        # write the solar file
    println("===== GML - Write Output File");
    name=string("results/Time", current_time, ".csv");
    cost = repeat([val_opt.cost], 12, 1)
    time = repeat([val_opt.time], 12, 1)
    global feeder_num=[string("feeder",1)]
    for i=2:12
        feeder_num=vcat(feeder_num, string("feeder",i))
    end
    RT_data_feeder=hcat(feeder_num, val_opt.Pf, val_opt.Qf, val_opt.Pg, val_opt.B,val_opt.R, cost, time)
    CSV.write(name, DataFrame(RT_data_feeder, [:Feeder, :Pf, :QF, :Pg, :B, :R, :Cost, :time]));
    # println("    ---- Finish writting files! ")
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
        RSRV_30=(ten_min_price), LMP_RT=(rt_lbmp), LMP_DA=(da_lbmp));
    return price_raw
end

# function read_demand_data()
#     filename = "../data/Feb_6_7.csv"
#     data_trace = CSV.read(filename)
#     pd_raw=zeros(12, 576)
#     for time=1:576
#         for feeder=1:12
#             pd_raw[feeder, time] = data_trace[time,feeder+1];
#         end
#     end
#     return pd_raw
# end
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
