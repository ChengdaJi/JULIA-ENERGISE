function GML_Sys_Ava(T, F, BN, pd, B_rate, R_rate, ancillary_type)

    println("GML - Start defining boundaries");
    ###############################################################################

    # feeder level
    Qf_max=0.05*positive_array(pd.traj);
    Qf_min = -Qf_max;
    # minimum solar
    Pg_min = zeros(F, T);
    # battery
    B_max = ones(F,1)*B_rate*ones(1,T)/12;
    B_min = zeros(F,T);
    R_max = R_rate*B_max;
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
    V_min = Base_V*0.96;
    V_max = Base_V*1.06;

    # impedan
    # repeat function
    r_c=[];
    x_c=[];
    for i=1:BN/3
        push!(r_c,r_1);
        push!(x_c,x_1);
    end
    r=vcat(r_c...);
    x=vcat(x_c...);
    # penalty
    beta=-15;
    println("GML - finish building boundaries")
    bd=bd_struct(T,F,BN,Pg_min,Qf_max,Qf_min,B_max,B_min,R_max,R_min,W,
    delta_t,P_rsrv_min,tau,k,V_max,V_min,r,x,beta,S);
    return bd
end

function optimal_stoach_scenario_offline(current_time, obj, pd, pg,
price, SN, B_feedback ,P_rsrv_feedback, ancillary_type);
    println("GML - start building constraint")

    #################
    #  F,T or BN, T #
    #################
    ## The parameters

    T = obj.T; # number of time slots
    F = obj.F; # number of feeders
    BN = obj.BN; # number of banks

    Pg_min = obj.Pg_min; # minimum power generation

    Qf_max = obj.Qf_max; # maximum reactive power generation
    Qf_min = obj.Qf_min; # minimum eactive power generation

    B_max = obj.B_max; # maximum storage level
    B_min = obj.B_min; # minimum storage level

    R_max = obj.R_max; # the maximum discharge rate
    R_min = obj.R_min; # the minimum discharge rate

    delta_t = obj.delta_t; # time interval

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
    TB=[];
    for tb_index=1:BN
        push!(TB, (tb_index-1)*4+1:tb_index*4)
    end
    V_min_0=69*0.96;
    V_max_0=69*1.04;
    icdf_095 = -1.6449;
    # A1_=[1 0 0; 0 1 0; 0 0 0];
    # c1_=[0; 0; 1];
    # A_=[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1];
    # c_=[0; 0; 1; 0]

    m = Model(with_optimizer(Mosek.Optimizer))
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

    println("===== Real Time Constraint Buildup =====")
    # @constraint(m, zeros(F,1).<=Pg_rt[1:F,1]);
    # @constraint(m, pg.mu_rt[:,1].>=Pg_rt[1:F,1]);
    # @constraint(m, R_min[:,1].<= R_rt[1:F,1]);
    # @constraint(m, R_rt[1:F,1].<= R_max[:,1]);
    # @constraint(m, B_rt[:,1].==B_feedback[:,current_time]);
    # @constraint(m, Qf_min[:,1].<=Qf_rt[1:F,1]);
    # @constraint(m, Qf_rt[1:F,1].<=Qf_max[:,1]);
    for feeder=1:F
        @constraint(m, 0<=Pg_rt[feeder,1]);
        @constraint(m, pg.mu_rt[feeder,1]>=Pg_rt[feeder,1]);
        @constraint(m, R_min[feeder,1]<= R_rt[feeder,1]);
        @constraint(m, R_rt[feeder,1]<= R_max[feeder,1]);
        @constraint(m, B_rt[feeder,1]==B_feedback[feeder,current_time]);
        @constraint(m, Qf_min[feeder,1]<=Qf_rt[feeder,1]);
        @constraint(m, Qf_rt[feeder,1]<=Qf_max[feeder,1]);
    end

    for bank=1:BN
        # temp_Pg_rt = Pg_rt[TB[bank][1],1];
        # temp_R_rt = R_rt[TB[bank][1],1];
        # temp_Qf_rt = Qf_rt[TB[bank][1],1];
        # for feeder_ite = 2:length(TB[bank])
        #     temp_Pg_rt = temp_Pg_rt+Pg_rt[TB[bank][feeder_ite],1];
        #     temp_R_rt = temp_R_rt+R_rt[TB[bank][feeder_ite],1];
        #     temp_Qf_rt = temp_Qf_rt+Qf_rt[TB[bank][feeder_ite],1];
        # end
        # @constraint(m, P_rt[bank,1]==
        #     sum(Pd[TB[bank][1]:TB[bank][end],1]
        #     -Pg_rt[TB[bank][1]:TB[bank][end],1]
        #     -R_rt[TB[bank][1]:TB[bank][end],1]));
        @constraint(m, P_rt[bank,1]==
            sum(Pd[TB[bank][1]:TB[bank][end],1])-
            sum(Pg_rt[TB[bank][feeder],1]-R_rt[TB[bank][feeder],1] for feeder=1:4));
        # @constraint(m, P_rt[bank,1]==
        #     sum(Pd[TB[bank][1]:TB[bank][end],1])
        #     -temp_Pg_rt-temp_R_rt);
        # @constraint(m, Q_rt[bank,1]==
        #     sum(Qf_rt[TB[bank][1]:TB[bank][end],1]));
        # @constraint(m, Q_rt[bank,1]==temp_Qf_rt);
        @constraint(m, Q_rt[bank,1]==sum(Qf_rt[TB[bank][feeder_ite],1] for feeder_ite=1:4));
        # @constraint(m, norm(A1_*[P_rt[bank,1]; Q_rt[bank,1]; S[1,1]])
        #     <=dot(c1_,[P_rt[bank,1];Q_rt[bank,1]; S[1,1]]));
        @constraint(m, [0.5*S[1,1], S[1,1], P_rt[bank,1], Q_rt[bank,1]] in RotatedSecondOrderCone());
        @constraint(m,P_hat_rt[bank,1]==
            r[bank]*l_rt[bank,1]+P_rt[bank,1]);
        @constraint(m,Q_hat_rt[bank,1]==
            x[bank]*l_rt[bank,1]+Q_rt[bank,1]);
        @constraint(m,v_rt[bank,1]==v_0_rt[1,1]+
            (r[bank]^2+x[bank]^2)*l_rt[bank,1]-
            2*(r[bank]*P_hat_rt[bank,1]+x[bank]*Q_hat_rt[bank,1]));
        # @constraint(m, norm(A_*[P_hat_rt[bank,1]; Q_hat_rt[bank,1];
        #     u_temp_rt[bank,1]; v_temp_rt[bank,1]])<=
        #     dot(c_,[P_hat_rt[bank,1]; Q_hat_rt[bank,1]; u_temp_rt[bank,1]; v_temp_rt[bank,1]]));
        @constraint(m, [0.5*l_rt[bank,1], v_0_rt[1,1], P_hat_rt[bank,1], Q_hat_rt[bank,1]] in RotatedSecondOrderCone());
        # @constraint(m, u_temp_rt[bank,1]==0.5*(v_0_rt[1,1]+l_rt[bank,1]))
        # @constraint(m, v_temp_rt[bank,1]==0.5*(v_0_rt[1,1]-l_rt[bank,1]))
        @constraint(m, V_min^2<= v_rt[bank,1]);
        @constraint(m, v_rt[bank,1]<= V_max^2);
    end
    # @constraint(m, V_min^2*ones(BN,1).<= v_rt[1:BN,1]);
    # @constraint(m, v_rt[1:BN,1].<= V_max^2*ones(BN,1));
    @constraint(m, v_0_rt[1,1]>=V_min^2);
    @constraint(m, v_0_rt[1,1]<=V_max^2);
    @constraint(m, sum(Q_hat_rt)<=0.05*sum(P_hat_rt));
    println("===== Stochastic Constraint Buildup =====")
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
    # @variable(m, u_temp[1:SN, 1:BN, 1:T-1]); # power distributed from bank
    # @variable(m, v_temp[1:SN, 1:BN, 1:T-1]); # reactive power distributed from bank

    for scenario = 1:SN
        # @constraint(m, Pg_min[:, 2:T].<=Pg[scenario, : ,:]);
        # @constraint(m, Pg[scenario, : ,:].<=
        #     positive_array(icdf_095.*sqrt.(pd.sigma+pg.sigma)+pg.mu_scenario));
        # @constraint(m, B_min[:,2:T] .<= B[scenario,:,:]);
        # @constraint(m, B[scenario,:,:].<= B_max[:,2:T]);
        # @constraint(m, R_min[:,2:T] .<= R[scenario,:,:]);
        # @constraint(m, R[scenario,:,:].<= R_max[:,2:T]);
        # for t=1:T-1
        #     if t==1
        #         @constraint(m, B[scenario,:,1] .== B_rt[:,1]-delta_t*R_rt[:,1])
        #     else
        #         @constraint(m, B[scenario,:,1] .==
        #             B[scenario,:,t-1] - R[scenario,:,t-1]*delta_t)
        #     end
        # end
        # @constraint(m, Qf_min[:,2:T] .<= Qf[scenario,:,:]);
        # @constraint(m, Qf[scenario,:,:].<= Qf_max[:,2:T]);
        for t=1:T-1
            for feeder=1:F
                @constraint(m, Pg_min[feeder, t]<=Pg[scenario, feeder ,t]);
                @constraint(m, Pg[scenario, feeder ,t]<=
                    positive_scalar(icdf_095*sqrt(pd.sigma[feeder,t]+pg.sigma[feeder,t])+pg.mu_scenario[feeder,t]));
                @constraint(m, B_min[feeder,t+1] <= B[scenario,feeder,t]);
                @constraint(m, B[scenario,feeder,t]<= B_max[feeder,t+1]);
                @constraint(m, R_min[feeder,t+1] <= R[scenario,feeder,t]);
                @constraint(m, R[scenario,feeder,t]<= R_max[feeder,t+1]);
                if t==1
                    @constraint(m, B[scenario,feeder,1] == B_rt[feeder,1]-delta_t*R_rt[feeder,1])
                else
                    @constraint(m, B[scenario,feeder,1] ==
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
                # @constraint(m, P_[scenario,bank,t]==
                #     sum(Pd[TB[bank][1]:TB[bank][end],t]
                #     -Pg[scenario,TB[bank][1]:TB[bank][end],t]
                #     -R[scenario,TB[bank][1]:TB[bank][end],t]));
                @constraint(m, P_[scenario,bank,t]==
                    sum(Pd[TB[bank][1]:TB[bank][end],t])
                    -temp_Pg-temp_R);
                # @constraint(m, Q_[scenario,bank,t]==
                #     sum(Qf[scenario,TB[bank][1]:TB[bank][end],t]));
                @constraint(m, Q_[scenario,bank,t]==temp_Qf);
                # @constraint(m, Q_[scenario,bank,t]==sum(Qf[scenario,TB[bank][feeder],t] for feeder=1:4));
                # @constraint(m, norm(A1_*[P_[scenario,bank,t]; Q_[scenario,bank,t]; S[1,t+1]])
                #     <=dot(c1_,[P_[scenario,bank,t];Q_[scenario,bank,t]; S[1,t+1]]));
                @constraint(m, [0.5*S[1,t+1], S[1,t+1], P_[scenario,bank,t], Q_[scenario,bank,t]] in RotatedSecondOrderCone());
                @constraint(m,P_hat[scenario,bank,t]==
                    r[bank]*l[scenario,bank,t]+P_[scenario,bank,t]);
                @constraint(m,Q_hat[scenario,bank,t]==
                    x[bank]*l[scenario,bank,t]+Q_[scenario,bank,t]);
                @constraint(m,v[scenario,bank,t]==v_0[scenario,t]+
                    (r[bank]^2+x[bank]^2)*l[scenario,bank,t]-
                    2*(r[bank]*P_hat[scenario,bank,t]+x[bank]*Q_hat[scenario,bank,t]));
                # @constraint(m, norm(A_*[P_hat[scenario,bank,t]; Q_hat[scenario,bank,t];
                #     u_temp[scenario,bank,t]; v_temp[scenario,bank,t]])<=
                #     dot(c_,[P_hat[scenario,bank,t]; Q_hat[scenario,bank,t];
                #      u_temp[scenario,bank,t]; v_temp[scenario,bank,t]]));
                @constraint(m, [0.5*l[scenario,bank,t], v_0[scenario,t], P_hat[scenario,bank,t], Q_hat[scenario,bank,t]] in RotatedSecondOrderCone());
                # @constraint(m, u_temp[scenario,bank,t]==0.5*(v_0[scenario,t]+l[scenario,bank,t]))
                # @constraint(m, v_temp[scenario,bank,t]==0.5*(v_0[scenario,t]-l[scenario,bank,t]))
                @constraint(m, V_min^2<=v[scenario,bank,t]);
                @constraint(m, V_max^2>=v[scenario,bank,t]);
            end
            @constraint(m, V_min_0^2<=v_0[scenario,t]);
            @constraint(m, V_max_0^2>=v_0[scenario,t]);
            @constraint(m, sum(Q_hat[scenario,:,t])<=0.05*sum(P_hat[scenario,:,t]));
        end
    end
    # @objective(m, Min, fn_cost_RHC(delta_t,price,P_hat,P_hat_rt,Pg,Pg_rt,beta,SN))
    @objective(m, Min,
        fn_cost_RHC(delta_t,price,P_hat_rt,P_hat,Pg_rt,Pg,pg,pd,beta,SN,icdf_095))
    status=optimize!(m);
    # println(status)
    println(termination_status(m))
    cost = JuMP.objective_value(m);
    time_solve=MOI.get(m, MOI.SolveTime());
    println(time_solve)
    println(cost)
end

function fn_cost_RHC(delta_t,price,P_hat_rt,P_hat,Pg_rt,Pg,pg,pd,beta,SN,icdf_095)
    println("load the optimal function")
    sum_prob = sum(price.probability[1:SN])
    P_hat_scenario=
        price.probability[1]/sum_prob*sum(P_hat[1, :, :], dims=1);
    Cost_P_hat_scenario = P_hat_scenario*reshape(price.lambda_scenario[1, :],287,1);
    Pg_diff_scenario=sum(Pg_rt[:,1])-sum(pg.mu_rt)+
        price.probability[1]/sum_prob*sum(Pg[1, :, :]);
    if SN>1
        for scenario = 2:SN
            P_hat_scenario=P_hat_scenario+
                price.probability[scenario]/sum_prob*sum(P_hat[scenario, :, :], dims=1);
            Cost_P_hat_scenario=Cost_P_hat_scenario+price.probability[scenario]/sum_prob*sum(P_hat[scenario, :, :], dims=1)*
                reshape(price.lambda_scenario[scenario, :],287,1);
            Pg_diff_scenario = Pg_diff_scenario+
                price.probability[scenario]/sum_prob*sum(Pg[scenario, :, :]);
        end
    end
    Pg_diff_scenario = Pg_diff_scenario-
        sum(positive_array(icdf_095.*sqrt.(pd.sigma+pg.sigma)+pg.mu_scenario));
    # println(length(Pg_diff_scenario))
    # P_hat_scenario=
    #     price.probability[1]*sum(P_hat[1, :, :], dims=1)+
    #     price.probability[2]*sum(P_hat[2, :, :], dims=1)+
    #     price.probability[3]*sum(P_hat[3, :, :], dims=1)+
    #     price.probability[4]*sum(P_hat[4, :, :], dims=1)+
    #     price.probability[5]*sum(P_hat[5, :, :], dims=1)+
    #     price.probability[6]*sum(P_hat[6, :, :], dims=1);
    # Pg_diff_scenario=sum(Pg_rt[:,1])-sum(pg.mu_rt)+
    #     price.probability[1]*sum(Pg[1, :, :])+
    #     price.probability[2]*sum(Pg[2, :, :])+
    #     price.probability[3]*sum(Pg[3, :, :])+
    #     price.probability[4]*sum(Pg[4, :, :])+
    #     price.probability[5]*sum(Pg[5, :, :])+
    #     price.probability[6]*sum(Pg[6, :, :])-
    #     sum(positive_array(icdf_095.*sqrt.(pd.sigma+pg.sigma)+pg.mu_scenario));
    # println(size(P_hat_scenario))
    # println(size(price.lambda_scenario[1, :]))
    # println(size(P_hat_scenario*reshape(price.lambda_scenario[1, :],287,1)))
    # println(size(Cost_P_hat_scenario))
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
