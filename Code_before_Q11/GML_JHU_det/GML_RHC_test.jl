function GML_Sys_Ava(T, F, BN, SN, pd, ancillary_type, icdf, B_cap)

    println("===== GML - Boundaries Buildup");

    ###############################################################################

    # feeder level
    Qf_max=0.05*positive_array(pd.traj);
    Qf_min = -Qf_max;
    # minimum solar
    Pg_min = zeros(F, T);
    # battery
    B_rate=B_cap;
    R_rate=1/3;
    B_max = ones(F,1)*B_rate*ones(1,T)/12;
    B_min = zeros(F,T);
    R_max = R_rate*B_max;
    R_min = -R_max;
    W=zeros(F,T);
    # ancillary

    tau=2;

    P_rsrv_min=zeros(1,T);
    k=12;

    delta_t = 1/12;
    S=35*ones(1,T);
    # V
    Base_V=69;
    V_min = Base_V*0.96;
    V_max = Base_V*1.06;

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

function optimal_stoach_scenario(current_time, obj, feedback, pd, pg, price, ancillary_type)
    println("===== GML - Optimization ")

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
    V_min_0=69*0.96;
    V_max_0=69*1.04;



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
        @constraint(m, Pg_rt[feeder,1]<=
            positive_scalar(icdf*sqrt(pd.sigma[feeder,1]+pg.sigma[feeder,1])+pg.mu[feeder,1]));
        @constraint(m, R_min[feeder,1]<= R_rt[feeder,1]);
        @constraint(m, R_rt[feeder,1]<= R_max[feeder,1]);
        @constraint(m, B_rt[feeder,1]==B_feedback[feeder,1]);
        @constraint(m, Qf_min[feeder,1]<=Qf_rt[feeder,1]);
        @constraint(m, Qf_rt[feeder,1]<=Qf_max[feeder,1]);
    end

    for bank=1:BN
        @constraint(m, P_rt[bank,1]==
            sum(Pd[TB[bank][1]:TB[bank][end],1])-
            sum(Pg_rt[TB[bank][feeder],1]+R_rt[TB[bank][feeder],1] for feeder=1:4));
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
    # @constraint(m, sum(Q_hat_rt)<=0.05*sum(P_hat_rt));



    # println(" ---- Stochastic Constraint Buildup")
    @variable(m, Pg[1:SN, 1:F, 1:T-1]); # the real power output
    @variable(m, Qf[1:SN, 1:F, 1:T-1]); # the real power output
    @variable(m, B[1:SN, 1:F, 1:T-1]); # the storage
    @variable(m, R[1:SN, 1:F, 1:T-1]);# the charge/discharge rate
    # bank level
    @variable(m, P_hat[1:SN, 1:BN, 1:T-1]);# the total delivered power on line
    @variable(m, Q_hat[1:SN, 1:BN, 1:T-1]);# the total delivered reactive on line
    @variable(m, l[1:SN, 1:BN, 1:T-1]);# current square
    @variable(m, v[1:SN, 1:BN, 1:T-1]);# voltage
    @variable(m, v_0[1:SN,1:T-1]);# voltage for P0
    @variable(m, P_[1:SN, 1:BN, 1:T-1]); # power distributed from bank
    @variable(m, Q_[1:SN, 1:BN, 1:T-1]); # reactive power distributed from bank

    for scenario = 1:SN
        for t=1:T-1
            for feeder=1:F
                @constraint(m, Pg_min[feeder, t]<=Pg[scenario, feeder ,t]);
                @constraint(m, Pg[scenario, feeder ,t]<=
                    positive_scalar(icdf*sqrt(pd.sigma[feeder,t+1]+pg.sigma[feeder,t+1])+pg.mu[feeder,t+1]));
                @constraint(m, B_min[feeder,t+1] <= B[scenario,feeder,t]);
                @constraint(m, B[scenario,feeder,t]<= B_max[feeder,t+1]);
                @constraint(m, R_min[feeder,t+1] <= R[scenario,feeder,t]);
                @constraint(m, R[scenario,feeder,t]<= R_max[feeder,t+1]);
                if t==1
                    @constraint(m, B[scenario,feeder,1] == B_rt[feeder,1]-delta_t*R_rt[feeder,1])
                else
                    @constraint(m, B[scenario,feeder,t] ==
                        B[scenario,feeder,t-1] - R[scenario,feeder,t-1]*delta_t)
                end
                @constraint(m, Qf_min[feeder,t+1] <= Qf[scenario,feeder,t]);
                @constraint(m, Qf[scenario,feeder,t]<= Qf_max[feeder,t+1]);
            end
            for bank=1:BN
                temp_Pg = Pg[scenario,TB[bank][1],t];
                temp_R = R[scenario,TB[bank][1],t];
                temp_Qf = Qf[scenario,TB[bank][1],t];
                for feeder_ite = 2:length(TB[bank])
                    temp_Pg = temp_Pg +Pg[scenario,TB[bank][feeder_ite],t];
                    temp_R = temp_R+R[scenario,TB[bank][feeder_ite],t];
                    temp_Qf = temp_Qf+Qf[scenario,TB[bank][feeder_ite],t];
                end
                @constraint(m, P_[scenario,bank,t]==
                    sum(Pd[TB[bank][1]:TB[bank][end],t])
                    -temp_Pg-temp_R);
                @constraint(m, Q_[scenario,bank,t]==temp_Qf);
                @constraint(m, [0.5*S[1,t+1], S[1,t+1], P_[scenario,bank,t], Q_[scenario,bank,t]] in RotatedSecondOrderCone());
                @constraint(m,P_hat[scenario,bank,t]==
                    r[bank]*l[scenario,bank,t]+P_[scenario,bank,t]);
                @constraint(m,Q_hat[scenario,bank,t]==
                    x[bank]*l[scenario,bank,t]+Q_[scenario,bank,t]);
                @constraint(m,v[scenario,bank,t]==v_0[scenario,t]+
                    (r[bank]^2+x[bank]^2)*l[scenario,bank,t]-
                    2*(r[bank]*P_hat[scenario,bank,t]+x[bank]*Q_hat[scenario,bank,t]));
                @constraint(m, [0.5*l[scenario,bank,t], v_0[scenario,t], P_hat[scenario,bank,t], Q_hat[scenario,bank,t]] in RotatedSecondOrderCone());
                @constraint(m, V_min^2<=v[scenario,bank,t]);
                @constraint(m, V_max^2>=v[scenario,bank,t]);
            end
            @constraint(m, V_min_0^2<=v_0[scenario,t]);
            @constraint(m, V_max_0^2>=v_0[scenario,t]);
            # @constraint(m, sum(Q_hat[scenario,:,t])<=0.05*sum(P_hat[scenario,:,t]));
        end
    end
    ###########################################################################
    # variables and Constraints for reserve markets
    if ancillary_type == "10min"
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
            temp_f_rsrv_c_fb = mult_fb[1]*P_rsrv_feedback[ini_fb]
            if length_fb >= 2
                for f_rsrv_fb_n=2:length_fb;
                    temp_f_rsrv_c_fb = temp_f_rsrv_c_fb+
                        mult_fb[f_rsrv_fb_n]*P_rsrv_feedback[ini_fb-1+f_rsrv_fb_n];
                end
            end
            @constraint(m, B_rsrv_rt==delta_t*temp_f_rsrv_c_fb)
        end
        @constraint(m, B_rsrv_rt <= sum(B_rt))
        # Scenario
        @variable(m, P_rsrv[1:SN,1:T-1])
        @variable(m, B_rsrv[1:SN,1:T-1])
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
                    @constraint(m, B_rsrv[scenario, t_real-current_time]==delta_t*temp_f_rsrv_c_fb);
                elseif fin == current_time
                    if current_time == 1
                        @constraint(m, B_rsrv[scenario, t_real-current_time]==delta_t*k*P_rsrv_rt);
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
                        @constraint(m, B_rsrv[scenario, t_real-current_time]==delta_t*k*P_rsrv_rt+delta_t*temp_f_rsrv_c_fb);
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
                                temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
                            end
                        end
                        @constraint(m, B_rsrv[scenario, t_real-current_time] == delta_t*mult_rt*P_rsrv_rt+delta_t*temp_f_rsrv_c_sc);
                    else
                        ini_sc = current_time+1;
                        fin_sc = fin;
                        length_sc = fin_sc - ini_sc +1;
                        mult_sc = k-length_sc+1:k;
                        temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time]
                        if length_sc > 1
                            for f_rsrv_sc_n=2:length_sc
                                temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
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
                        delta_t*mult_rt*P_rsrv_rt+delta_t*temp_f_rsrv_c_sc+delta_t*temp_f_rsrv_c_fb);
                    end
                elseif ini == current_time
                    ini_sc = current_time+1;
                    fin_sc = fin;
                    mult_sc = 2:k;
                    temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time];
                    for f_rsrv_sc_n=2:k-1
                        temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
                    end
                    @constraint(m, B_rsrv[scenario, t_real-current_time] ==delta_t*P_rsrv_rt+delta_t*temp_f_rsrv_c_sc);
                elseif ini > current_time
                    ini_sc = ini;
                    fin_sc = fin;
                    mult_sc = 1:k;
                    temp_f_rsrv_c_sc = mult_sc[1]*P_rsrv[scenario, ini_sc-current_time]
                    for f_rsrv_sc_n=2:k
                        temp_f_rsrv_c_sc=temp_f_rsrv_c_sc+mult_sc[f_rsrv_sc_n]*P_rsrv[scenario, ini_sc-1+f_rsrv_sc_n-current_time];
                    end
                    @constraint(m, B_rsrv[scenario, t_real-current_time] == delta_t*temp_f_rsrv_c_sc);
                    if  t_real-current_time >= T-tau-1
                        @constraint(m, P_rsrv[scenario, t_real-current_time]==0)
                    end
                end
                @constraint(m, B_rsrv[scenario, t_real-current_time] <= sum(B[scenario, :, t_real-current_time]))
            end
        end
    end
    if ancillary_type == "10min"
        @objective(m, Min,
            fn_cost_RHC_anc(delta_t,P_hat_rt,P_hat,Pg_rt,Pg,P_rsrv_rt,P_rsrv,price,pg,pd,beta,SN,obj))
    elseif ancillary_type == "without"
        @objective(m, Min,
            fn_cost_RHC_rt(delta_t,P_hat_rt,P_hat,Pg_rt,Pg,price,pg,pd,beta,SN,obj))
    else
        println("Invalid ancillary type");
        exit()
    end
    status=optimize!(m);
    cost_printout = JuMP.objective_value(m);
    time_solve=MOI.get(m, MOI.SolveTime());
    println(string("    ----Solve Time: ", time_solve))
    println(string("    ----Optimal Cost for whole horion: ", cost_printout))
    ## obtaining value
    Qf_o=JuMP.value.(Qf_rt)
    Pg_o=JuMP.value.(Pg_rt)
    B_o=JuMP.value.(B_rt)
    R_o=JuMP.value.(R_rt)
    P_hat_o=JuMP.value.(P_hat_rt)
    Q_hat_o=JuMP.value.(Q_hat_rt)
    l_o=JuMP.value.(l_rt)
    v_o=JuMP.value.(v_rt)
    v_0_o=JuMP.value.(v_0_rt)
    if ancillary_type == "10min" || ancillary_type == "30min"
        P_rsrv_o=JuMP.value(P_rsrv_rt);
        B_rsrv_o=JuMP.value(B_rsrv_rt);
        P_rsrv_s=JuMP.value.(P_rsrv);
        B_rsrv_s=JuMP.value.(B_rsrv);
    else
        P_rsrv_o = 0;
        B_rsrv_o = 0;
    end
    # in S 1
    Pg_s = JuMP.value.(Pg[1,:,:])

    P_0_o = calculate_P0_real(pg, pd, price, R_o, l_o, obj, TB);
    # println(string("old p0", sum(P_hat_o)))
    # println(string("new p0", P_0_o))
    if ancillary_type == "without"
        cost_o=P_0_o*price.lambda_ct/12 + beta*(sum(Pg_o)-sum(pg.mu_ct))/12;
    elseif ancillary_type == "10min" || ancillary_type == "30min"
        cost_o=P_0_o*price.lambda_ct/12 + beta*(sum(Pg_o)-sum(pg.mu_ct))/12-
        delta_t*price.alpha_ct*P_rsrv_o;
    end
    println(string("    ----Optimal Cost at current time: ", cost_o))
    Pf_o=zeros(F,1)
    for feeder = 1:F
        Pf_o[feeder,1]=Pd[feeder,1]-Pg_o[feeder,1]-R_o[feeder,1]
    end

    val_opt = (Pf = (Pf_o), Qf=(Qf_o), Pg = (Pg_o), B=(B_o), R=(R_o), P_hat = (P_hat_o), Q_hat = (Q_hat_o),
    l = (l_o), v = (v_o), v_0=(v_0_o), P_rsrv = (P_rsrv_o), B_rsrv = (B_rsrv_o), P_0 = (P_0_o),
    cost=(cost_o), time = (time_solve))
    return val_opt
end

function fn_cost_RHC_anc(delta_t,P_hat_rt,P_hat,Pg_rt,Pg,P_rsrv_rt,P_rsrv,price,pg,pd,beta,SN,obj)

    # println("    ---- Load the optimal function")
    T=obj.T;
    icdf = obj.icdf;
    sum_prob = sum(price.probability[1:SN])


    lambda_ct=price.probability[1]/sum_prob*price.lambda_scenario[1,1];
    alpha_ct=price.probability[1]/sum_prob*price.alpha_scenario[1,1]
    if SN>1
        for scenario = 2:SN
            lambda_ct=lambda_ct+
                price.probability[scenario]/sum_prob*price.lambda_scenario[scenario,1];
            alpha_ct=alpha_ct+
                price.probability[scenario]/sum_prob*price.alpha_scenario[scenario,1];
            end
    end
    # Current time
    Cost_P_hat_ct = delta_t*lambda_ct*sum(P_hat_rt[:,1]);

    Pg_diff_ct = sum(Pg_rt[:,1]) -
        sum(positive_array(icdf.*sqrt.(pd.sigma[:,1]+pg.sigma[:,1])+pg.mu[:,1]));
    Cost_Pg_diff_ct = delta_t*beta*Pg_diff_ct;

    Revenue_P_rsrv_ct = delta_t*alpha_ct*P_rsrv_rt;

    # Future
    P_hat_scenario=
        price.probability[1]/sum_prob*sum(P_hat[1, :, :], dims=1);
    Cost_P_hat_scenario = delta_t*P_hat_scenario*reshape(price.lambda_scenario[1, 2:end],T-1,1);

    Pg_scenario=price.probability[1]/sum_prob*sum(Pg[1, :, :]);

    P_rsrv_scenario = price.probability[1]/sum_prob*P_rsrv[1,:];
    Revenue_P_rsrv_scenario = delta_t*reshape(P_rsrv_scenario, 1, T-1)*
        reshape(price.alpha_scenario[1, 2:end],T-1,1);
    if SN>1
        for scenario = 2:SN
            P_hat_scenario=P_hat_scenario+
                price.probability[scenario]/sum_prob*sum(P_hat[scenario, :, :], dims=1);
            Cost_P_hat_scenario=Cost_P_hat_scenario+delta_t*price.probability[scenario]/sum_prob*sum(P_hat[scenario, :, :], dims=1)*
                reshape(price.lambda_scenario[scenario, 2:end],T-1,1);
            Pg_scenario = Pg_scenario+
                price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :]);
            P_rsrv_scenario = price.probability[scenario]/sum_prob*P_rsrv[scenario,:];
            Revenue_P_rsrv_scenario = Revenue_P_rsrv_scenario+delta_t*
                reshape(P_rsrv_scenario, 1, T-1)*reshape(price.alpha_scenario[scenario, 2:end],T-1,1);
        end
    end

    # FOL cost part
    Pg_diff_scenario = Pg_scenario-
        sum(positive_array(icdf.*sqrt.(pd.sigma[:,2:end]+pg.sigma[:,2:end])
        +pg.mu[:,2:end]));
    Cost_Pg_diff_scenario = delta_t*beta*Pg_diff_scenario;
    println(string("    ----Case: Real-time Balancing and 10 min Reserve Market"))
    return Cost_P_hat_ct + Cost_Pg_diff_ct- Revenue_P_rsrv_ct +
            Cost_P_hat_scenario[1,1] + Cost_Pg_diff_scenario - Revenue_P_rsrv_scenario[1,1]
end

function fn_cost_RHC_rt(delta_t,P_hat_rt,P_hat,Pg_rt,Pg,price,pg,pd,beta,SN,obj)

    # println("    ---- Load the optimal function")
    T=obj.T;
    icdf = obj.icdf;
    sum_prob = sum(price.probability[1:SN])
    P_hat_scenario=
        price.probability[1]/sum_prob*sum(P_hat[1, :, :], dims=1);
    Cost_P_hat_scenario = P_hat_scenario*reshape(price.lambda_scenario[1, :],T-1,1);
    Pg_diff_scenario=sum(Pg_rt[:,1])-sum(pg.mu_rt)+
        price.probability[1]/sum_prob*sum(Pg[1, :, :]);
    if SN>1
        for scenario = 2:SN
            P_hat_scenario=P_hat_scenario+
                price.probability[scenario]/sum_prob*sum(P_hat[scenario, :, :], dims=1);
            Cost_P_hat_scenario=Cost_P_hat_scenario+price.probability[scenario]/sum_prob*sum(P_hat[scenario, :, :], dims=1)*
                reshape(price.lambda_scenario[scenario, :],T-1,1);
            Pg_diff_scenario = Pg_diff_scenario+
                price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :]);
        end
    end
    Pg_diff_scenario = Pg_diff_scenario-
        sum(positive_array(icdf.*sqrt.(pd.sigma+pg.sigma)+pg.mu_scenario));
    println("    ----Case: Real-time Balancing")
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

function write_output_out(val_opt, filename)
        # write the solar file
    println("===== GML - Write Output File");
    # name=string("results/Solar0025Time", current_time, ".csv");
    cost = repeat([val_opt.cost], 12, 1)
    time = repeat([val_opt.time], 12, 1)
    P_0 = repeat([val_opt.P_0], 12, 1)
    p_rsrv = repeat([val_opt.P_rsrv], 12, 1)
    B_rsrv = repeat([val_opt.B_rsrv], 12, 1)
    v_0 = zeros(12,1)
    v_0[1]=val_opt.v_0[1,1]
    v=zeros(12,1)
    v[1:3]=val_opt.v
    # v = [val_opt.v, zeros(9,1)]
    l=zeros(12,1)
    l[1:3]=val_opt.l;
    l = [val_opt.l; zeros(9,1)]
    P_hat=zeros(12,1)
    P_hat[1:3]=val_opt.P_hat;
    Q_hat=zeros(12,1)
    Q_hat[1:3]=val_opt.Q_hat;
    global feeder_num=[string("feeder",1)]
    for i=2:12
        feeder_num=vcat(feeder_num, string("feeder",i))
    end
    RT_data_feeder=hcat(feeder_num, val_opt.Pf, val_opt.Qf, val_opt.Pg, val_opt.B, val_opt.R, p_rsrv, B_rsrv, P_0, cost, time, v_0, v, l, P_hat, Q_hat)
    CSV.write(filename, DataFrame(RT_data_feeder, [:Feeder, :Pf, :QF, :Pg, :B, :R, :P_rsrv, :B_rsrv, :P_0, :Cost, :time, :v_0, :v, :l, :P_hat, :Q_hat]));
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

function calculate_P0_real(pg, pd, price, R, l, obj, TB)
    # calculate
    P_bank=zeros(obj.BN,1);
    P_line=zeros(obj.BN,1);
    for bank=1:obj.BN
        P_bank[bank,1]=
            sum(pd.ct[TB[bank][feeder],1]-pg.mu_ct[TB[bank][feeder],1]-R[TB[bank][feeder],1] for feeder=1:4);
        P_line[bank,1]=
            obj.r[bank]*l[bank,1]+P_bank[bank,1];
    end
    P_0 = sum(P_line[bank,1] for bank=1:obj.BN)
    return P_0
end
