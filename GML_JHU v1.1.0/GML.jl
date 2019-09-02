function GML(Pd, Sg_max, Lambda, Beta, Alpha)
    # Input:
    # -Pd: real-time demand data. Each row represents a feeder, whereas and each column represents a time slot,
    # e.g., in the small test case, there are 12 feeders with 288 slots (every 5 minutes for 24 hours). Therefore Pd is a 12 by 288 matrix.
    # -Pg_trend: solar data.  A vector. E.g., in our test case we have 2 location and 288 slots, therefore, Pg_trend a 2 by 288 matrix.
    # -Lambda: Real-time market price data. A vector, same form as Pg_trend.
    # -Beta: Solar (curtailment) penalty. A positive scalar.
    # -alpha: Ancillary market price data. A column vector, same form as Pg_trend.
    #
    # Output:
    # A solar CSV file "solar.csv" comprising the optimal solar trajectory for the FOL.
    #  	-column 1: timestamp.
    # 	-column 2 to column 13: optimal solar generation P_g at each feeder.
    # 	-column 14 to column 25: optimal reactive solar generation Q_g
    #
    # A battery CSV file "battery.csv" comprising the optimal battery trajectory for the FOL.
    #     -column 1: timestamp.
    #     -column 2 to column 13: power discharge R_f at each feeder
    #     -column 14 to column 25: reactive power discharge C_f at each feeder
    #     -column 25 to column 36: battery SOC B_f at each feeder


    println("GML - Start Optimization");

    ###############################################################################
    ### Generate boundaries
    bd=GML_boundaries(Pd, Sg_max, Lambda, Beta, Alpha);
    ###############################################################################
    ### optimization
    Opt_value=GML_optmization(bd);
    GML_output(Opt_value);
    return Opt_value.Pf
end

function GML_boundaries(Pd, Sg_max, Lambda, Beta, Alpha)
    # Function GML_boundaries(Pd, Pg_trend, Lambda, Beta, Alpha) is the function
     # that returns the required upper and lower bounds on .......UPDATE this
    # Inputs:
    # Same as GML();
    # Outputs:
    # -bd: composed by constructor "bd_struct" and it contains the following information:
    	#number of time slots T;
            #number of feeders F;
            #number of banks BN;
            # demand Pd;
            # reactive demand Qd;
            # maximum power output Pg_max;
            # minimum power output Pg_min;
            # maximum reactive power output Qg_max;
            # minimum reactive power output Qg_min;
            # the maximum power of storage B_max;
            # the minimum power of storage B_min;
            # the maximum power discharge rate R_max;
            # the minimum power discharge rate R_min;
            # power storage C_max;
            # the minimum power of storage C_min;
            # power storage b_0;
            # the exogenous change  W;
            # time interval delta_t;
            # The maximum ancillary power P_rsrv_min;
            # the apparent power S;
            # the maximum voltage V_max;
            # the minimum voltage V_min;
            # price of ancillary power alpha;
            # price of cost beta;
            # the grid spot price lambda;
            # ancillary respond time tau;
            # impedance r;
            # impedance  x;
            # time asso/w ancellary k;
    # Assumption on bounds:
    #  	Pg_max = 0.5 max(Pd)*Pg_trends/sum(Pg_trends)*sum(Pd)
    # 	Qg_max = 0.05*Pg_max
    # 	Pg_min = 0;
    # 	Qg_min = -Qg_max;
    # 	B_max = 0.1*Pd
    # 	B_min=0;
    # 	R_max = B_max;
    # 	R_min = - B_max;
    # 	C_max = 0.05*R_max;
    # 	C_min = -0.05*R_max;
    # 	B_f(0)=0.5*B_max_f(0);
    # 	B_f(t_f)=0.5*B_max_f(t_f);	(t_f is the last slot)
    # 	if B_max(t)-B_max(t-1)>0
    # 		W_f_t=0.1*(B_max(t)-B_max(t-1))
    # 	otherwise
    # 		W_f_t=0.9*(B_max(t)-B_max(t-1))
    # 	P_rsrv_min=0;

    println("GML - Start defining boundaries");
    ###############################################################################

    number_=size(Pd);
    F=number_[1];
    T=number_[2];
    Qd=0.05*Pd;

    maxPd=zeros(1,F)
    # maxQd=zeros(1,F)
    for i=1:F
        maxPd[i]=maximum(Pd[i,:])
        # maxQd[i]=maximum(Qd[i,:])
    end
    # Total_Pd1=sum(sum(Pd[1:8,:]));
    # Total_Pd2=sum(sum(Pd[9:F,:]));
    Sg_max = Sg_max
    ###############################################################################
    # generation

    # Pg_trends1=Pg_trends[:,1]';
    # Pg_trends2=Pg_trends[:,2]';
    # Total_Pg1=sum(sum(maxPd[1:8]*Pg_trends1));
    # Total_Pg2=sum(sum(maxPd[9:F]*Pg_trends2));


    Pg_min = zeros(F,T);
    # Pg_max = zeros(F,T);
    # Pg_max[1:8,:]=0.5*maxPd[1:8]*Pg_trends1/Total_Pg1*Total_Pd1;
    # Pg_max[9:F,:]=0.5*maxPd[9:F]*Pg_trends2/Total_Pg2*Total_Pd2;

    # Qg_max=0.05*Pg_max;
    # Qg_min=-Qg_max;
    ###############################################################################
    # storage

    B_max = 0.1*Pd;
    B_min = zeros(F,T);
    R_max = B_max;
    R_min = -B_max;
    C_max = 0.05*R_max;
    C_min = -0.05*R_max
    b_0 = 0.5*B_max[:,1]
    W=zeros(F,T);

    for i=1:F
        for t=2:T
            if B_max[i,t]-B_max[i,t-1]>=0
                W[i,t]=0.1*(B_max[i,t]-B_max[i,t-1]);
            else
                W[i,t]=0.9*(B_max[i,t]-B_max[i,t-1]);
            end
        end
    end

    ###############################################################################
    delta_t = 1/12;
    ###############################################################################
    P_rsrv_min=zeros(1,T)
    tau=2;
    ###############################################################################
    S=35*ones(1,T)
    Base_V=69;
    V_min = Base_V*0.96;
    V_max = Base_V*1.06;
    ## impedance
    r=[0.13275, 0.13275, 0.199125, 0.13275, 0.13275, 0.199125, 0.13275, 0.13275, 0.199125];
    x=[1.00426, 1.00426, 1.50639, 0.13275, 0.13275, 0.199125, 0.13275, 0.13275, 0.199125];
    k=12;
    BN=9;
    ###############################################################################
    bd_return=bd_struct(T, F, BN, Pd, Qd, Sg_max, Pg_min, B_max,
            B_min, R_max, R_min, C_max, C_min, b_0, W, delta_t,
            P_rsrv_min,S,V_max, V_min, Alpha, Beta, Lambda, tau, r, x, k);
    println("GML - Finish defining boundaries")
    return bd_return

end

function GML_optmization(bd)
    println("GML - Start performing optimization");
    #################
    #  F,T or BN, T #
    #################
    ## The parameters

    T = bd.T; # number of time slots
    F = bd.F; # number of feeders
    BN = bd.BN; # number of banks
    Pd = bd.Pd; # demand data (input)
    Qd = bd.Qd; # demand (input)
    Sg_max=bd.Sg_max;
    # Pg_max = bd.Pg_max; # maximum power generation
    Pg_min = bd.Pg_min; # minimum power generation

    # Qg_max = bd.Qg_max; # maximum reactive power generation
    # Qg_min = bd.Qg_min; # minimum eactive power generation

    B_max = bd.B_max; # maximum storage level
    B_min = bd.B_min; # minimum storage level

    R_max = bd.R_max; # the maximum discharge rate
    R_min = bd.R_min; # the minimum discharge rate

    C_max = bd.C_max; # maximum storage level
    C_min = bd.C_min; # minimum storage level

    b_0 = bd.b_0; # storage level
    W = bd.W; # exogenour change
    delta_t = bd.delta_t; # time interval

    P_rsrv_min = bd.P_rsrv_min; # The minimum ancillary power

    S=bd.S; # the apparents power

    V_min=bd.V_min;
    V_max=bd.V_max;
    alpha=bd.alpha; # price in ancillary power
    beta=bd.beta; # price in cost function
    lambda= bd.lambda; # the grid spot price
    tau = bd.tau;

    r=bd.r; # the resistance
    x=bd.x; # the reactance
    k=bd.k; # time with ancellary

    ###############################################################################
    # define variable
    m = Model(solver=GurobiSolver(Presolve=0))
    @variable(m, Pg[1:F, 1:T])
    @variable(m, Qg[1:F, 1:T])
    @variable(m, B[1:F, 1:T])
    @variable(m, R[1:F, 1:T])
    @variable(m, C[1:F, 1:T])
    @variable(m, P_hat[1:BN, 1:T])
    @variable(m, Q_hat[1:BN, 1:T])
    @variable(m, l[1:BN, 1:T])
    @variable(m, v[1:BN, 1:T])
    @variable(m, v_0[1, 1:T])
    @variable(m, P_[1:BN, 1:T])
    @variable(m, Q_[1:BN, 1:T])
    @variable(m, P_rsrv[1, 1:T])
    @variable(m, B_rsrv[1, 1:T])
    @variable(m, u_[1:BN, 1:T])
    @variable(m, v_[1:BN, 1:T])
    ###############################################################################
    # constraint 5a)-5d) on ACC paper
    @constraint(m, Pg.>= zeros(F,T))
    @constraint(m, Pg.>= Pg_min)
    @constraint(m, B.<= B_max)
    @constraint(m, B.>= B_min)
    @constraint(m, R.<= R_max)
    @constraint(m, R.>= R_min)
    @constraint(m, C.<= C_max)
    @constraint(m, C.>= C_min)
    ##############################################################################
    # constraint 5f) on ACC paper
    println("in constraint 5f")
    for t=1:T-1
        if t==1
            @constraint(m, B[:,t] .== b_0[:,t])
        end
        @constraint(m, B[:,t+1] .== B[:,t] - R[:,t]*delta_t + W[:,t])
    end
    @constraint(m, B[:,T] .== 0.5*B_max[:,T])
    ###############################################################################
    # constraint 5g)-5i) on ACC paper
    println("in constraint 5g-i")
    for t=1:T
        if t<=6
        @constraint(m, B_rsrv[1,t]==0)
        elseif t>T-tau+1
            start=max(t-k+1-tau,1);
            K=1:t-tau-start+1;
            # @constraint(m,P_rsrv[1,t]==0);
            @constraint(m,B_rsrv[1,t]==delta_t*K'*P_rsrv[1,start:t-tau]);
        else
           start=max(t-k+1-tau,1);
            K=1:t-tau-start+1;
            @constraint(m,B_rsrv[1,t]==delta_t*K'*P_rsrv[1,start:t-tau]);
        end
        @constraint(m, sum(B[:,t])>=B_rsrv[1,t])
    end
    @constraint(m, B_rsrv[1,T]==0)
    @constraint(m,P_rsrv[1,:] .>= P_rsrv_min[1,:])
    ###############################################################################
    # constraint 5i)-5j) on ACC paper
    println("in constraint 5i-j")
    TBl=[0 4 8 12 16 20 24 28 32 36]
    for bank=2:BN+1
        for t=1:T
             @constraint(m, P_[bank-1,t] == sum(Pd[TBl[bank-1]+1:TBl[bank],t]
                -Pg[TBl[bank-1]+1:TBl[bank],t]
                -R[TBl[bank-1]+1:TBl[bank],t]));
        end
    end

    for bank=2:BN+1
        for t=1:T
            @constraint(m, Q_[bank-1,t] == sum(Qd[TBl[bank-1]+1:TBl[bank],t]
                -Qg[TBl[bank-1]+1:TBl[bank],t]
                -C[TBl[bank-1]+1:TBl[bank],t]));
        end
    end

    ###############################################################################
    # constraint 5l) on ACC paper
    println("in constraint 5l")
    for bank=1:BN
        for t=1:T
            A1_=[1 0 0; 0 1 0; 0 0 0];
            c1_=[0; 0; 1];
            @constraint(m, norm(A1_*[P_[bank,t];Q_[bank,t]; S[1,t]])
                <=dot(c1_,[P_[bank,t];Q_[bank,t]; S[1,t]]));
        end
    end
    ###############################################################################
    # constraint 5m) on ACC paper
    println("in constraint 5m")
    @constraint(m, V_min^2*ones(BN,T) .<= v);
    @constraint(m, v .<= V_max^2*ones(BN,T));

    ############################
    V_min_0=69*0.96;
    V_max_0=69*1.04;
    #############################
    @constraint(m, V_min_0^2 .<= v_0[1,:]);
    @constraint(m, v_0[1,:] .<= V_max_0^2);

    ################################################################################
    # constraint 5n)-5o) on ACC paper
    println("in constraint 5n-0")
    for i=1:BN
        for t=1:T
            @constraint(m, P_hat[i,t]==r[i]*l[i,t]+P_[i,t]);
            @constraint(m,Q_hat[i,t]==x[i]*l[i,t]+Q_[i,t]);
        end
    end

    ################################################################################
    # constraint 5p)-5q) on ACC paper
    println("in constraint 5p-q")
    for i=1:BN
        for t=1:T
            @constraint(m, v[i,t]==v_0[1,t]+(r[i]^2+x[i]^2)*l[i,t]-
                2*(r[i]*P_hat[i,t]+x[i]*Q_hat[i,t]));
        end
    end

    # #
    for i=1:BN
        for t=1:T
            A_=[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1];
            c_=[0; 0; 1; 0]
            @constraint(m, norm(A_*[P_hat[i,t]; Q_hat[i,t]; u_[i,t]; v_[i,t]])<=
                dot(c_,[P_hat[i,t]; Q_hat[i,t]; u_[i,t]; v_[i,t]]));
        end
    end

    for i=1:BN
        for t=1:T
            @constraint(m, u_[i,t]==0.5*(v_0[1,t]+l[i,t]))
            @constraint(m, v_[i,t]==0.5*(v_0[1,t]-l[i,t]))
        end
    end

    ################################################################################
    # New Constraints on Solar
    println("in constraint solar")
    for t=1:T
        for f=1:F
        # println( size(Pg[:,t]));
        # println( size(sum(Pg[:,t])));
        @constraint(m, norm([Pg[f,t]; Qg[f,t]]) <= Sg_max[f,t]);
        end
    end
    ################################################################################

    @objective(m, Min, fn_cost(delta_t,lambda, P_hat, beta, Pg, Qg, Sg_max ,alpha,
    P_rsrv))
    status=solve(m);


    Pg_o=getvalue(Pg)
    Qg_o=getvalue(Qg)
    B_o=getvalue(B)
    R_o=getvalue(R)
    C_o=getvalue(C)
    P_hat_o=getvalue(P_hat)
    Q_hat_o=getvalue(Q_hat)
    l_o=getvalue(l)
    v_o=getvalue(v)
    P_rsrv_o=getvalue(P_rsrv)
    B_rsrv_o=getvalue(P_rsrv)
    P0_0=sum(P_hat_o)
    cost_o=getobjectivevalue(m);
    Pf_o=Pd-Pg_o-R_o;
    Qf_o=Qd-Qg_o-C_o;

    optimized_value=Optimization_output_struct(
    Pf_o,
    Qf_o,
    Pg_o,
    Qg_o,
    B_o,
    R_o,
    C_o,
    P_hat_o,
    Q_hat_o,
    l_o,
    v_o,
    P_rsrv_o,
    B_rsrv_o,
    P0_0,
    cost_o)
    println("GML - Finish performing optimization");
    println("GML - Optimal Cost is: ", cost_o, " USD");
    return optimized_value
end

###############################################################################
function fn_cost(delta_t,lambda,P_hat,beta,Pg, Qg, Sg_max ,alpha,P_rsrv)
    println("in optimal function")
    # This function define the operating cost
    println(size(sum(Pg)))
    PP=sum(P_hat,dims=1);
    println(size(delta_t*(lambda[:,1]'*PP'-alpha[:,:]'*P_rsrv[:])))
    # PP2=sum((Pg-Pg_max));
    return delta_t*(lambda[:,1]'*PP'
    -beta*(sum(Sg_max)-sum(Pg))
    -alpha[:,:]'*P_rsrv[:])[1,1]
end

###############################################################################
function GML_output(opt_value)
    # write the solar file
    println("GML - Start writing files!");
    total_n=size(opt_value.Pg)[2];
    time_stamp=[1:total_n]
    time_stamp=zeros(1,total_n)
    for i=1:size(opt_value.Pf)[2]
        time_stamp[i]=i
    end

    Power=[time_stamp; opt_value.Pf]';
    # CSV.write("Power.csv", DataFrame(Power, [:timestamp, :bnk139_Pf_1,
    # :bnk139_Pf_2, :bnk139_Pf_3, :bnk139_Pf_4, :bnk239_Pf_5, :bnk239_Pf_6,
    # :bnk239_Pf_7, :bnk239_Pf_8, :bnk276_Pf_9, :bnk276_Pf_10, :bnk276_Pf_11,
    # :bnk276_Pf_12]));
    # println("GML - Finish writting files! ")
end
