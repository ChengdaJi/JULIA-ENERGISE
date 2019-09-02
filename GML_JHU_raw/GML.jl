function GML(Pd, Pg_trend, Lambda, Beta, Alpha)
    println("GML - Start Optimization");

###############################################################################
### Generate boundaries
    bd=GML_boundaries(Pd, Pg_trend, Lambda, Beta, Alpha);
###############################################################################
### optimization
    Opt_value=GML_optmization(bd);
    GML_output(Opt_value);
    println("GML - Finish Optimization");
    println(Opt_value.cost)

end

function GML_boundaries(Pd, Pg_trends, Lambda, Beta, Alpha)
    println("GML - Start defining boundaries");
###############################################################################
    # demand
    # Pd=demand.Pd/1000
    number_=size(Pd)
    F=number_[1];
    T=number_[2];
    Qd=0.05*Pd;

    # Qd=demand.Qd/1000
    maxPd=zeros(1,F)
    # maxQd=zeros(1,F)
    for i=1:F
        maxPd[i]=maximum(Pd[i,:])
        # maxQd[i]=maximum(Qd[i,:])
    end
    Total_Pd1=sum(sum(Pd[1:8,:]));
    Total_Pd2=sum(sum(Pd[9:F,:]));
    # Total_Qd1=sum(sum(Qd[1:8,:]));
    # Total_Qd2=sum(sum(Qd[9:F,:]));
###############################################################################
    # generation

    Pg_trends1=Pg_trends[:,1]';
    Pg_trends2=Pg_trends[:,2]';
    # println(size(Pg_trends1))
    # println(size(maxPd[1:8]))
    # println(maxPd[1:8]*Pg_trends1)
    Total_Pg1=sum(sum(maxPd[1:8]*Pg_trends1));
    Total_Pg2=sum(sum(maxPd[9:F]*Pg_trends2));


    Pg_min = zeros(F,T);
    Pg_max = zeros(F,T);
    Pg_max[1:8,:]=0.5*maxPd[1:8]*Pg_trends1/Total_Pg1*Total_Pd1;
    Pg_max[9:F,:]=0.5*maxPd[9:F]*Pg_trends2/Total_Pg2*Total_Pd2;

    Qg_max=0.05*Pg_max;
    Qg_min=-Qg_max;
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
    r=[0.13275, 0.13275, 0.199125];
    x=[1.00426, 1.00426, 1.50639];
    k=12;
###############################################################################



    bd_return=bd_struct(T, F, 3, Pd, Qd, Pg_max, Pg_min, Qg_max, Qg_min, B_max,
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

    Pg_max = bd.Pg_max; # maximum power generation
    Pg_min = bd.Pg_min; # minimum power generation

    Qg_max = bd.Qg_max; # maximum reactive power generation
    Qg_min = bd.Qg_min; # minimum eactive power generation

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
    @constraint(m, Pg.>= zeros(F,T))
    @constraint(m, Pg.<= Pg_max)
    @constraint(m, Pg.>= Pg_min)
    @constraint(m, Qg.<= Qg_max)
    @constraint(m, Qg.>= Qg_min)
    @constraint(m, B.<= B_max)
    @constraint(m, B.>= B_min)
    @constraint(m, R.<= R_max)
    @constraint(m, R.>= R_min)
    @constraint(m, C.<= C_max)
    @constraint(m, C.>= C_min)
# #
# # ##############################################################################
# #
    for t=1:T-1
        if t==1
            @constraint(m, B[:,t] .== b_0[:,t])
        end
        @constraint(m, B[:,t+1] .== B[:,t] - R[:,t]*delta_t + W[:,t])
    end
    @constraint(m, B[:,T] .== 0.5*B_max[:,T])
# #
# # ###############################################################################
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
# # ###############################################################################
    @constraint(m,P_rsrv[1,:] .>= P_rsrv_min[1,:])
# # ###############################################################################
TBl=[0 4 8 12]
#
# #
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

# # ###############################################################################
# #
for bank=1:BN
    for t=1:T
        A1_=[1 0 0; 0 1 0; 0 0 0];
        c1_=[0; 0; 1];
        @constraint(m, norm(A1_*[P_[bank,t];Q_[bank,t]; S[1,t]])
            <=dot(c1_,[P_[bank,t];Q_[bank,t]; S[1,t]]));
    end
end

# #
@constraint(m, V_min^2*ones(BN,T) .<= v);
@constraint(m, v .<= V_max^2*ones(BN,T));

# # # ############################
V_min_0=69*0.96;
V_max_0=69*1.04;
# # #############################
@constraint(m, V_min_0^2 .<= v_0[1,:]);
@constraint(m, v_0[1,:] .<= V_max_0^2);

# #
# #
# # ################################################################################
for i=1:BN
    for t=1:T
        @constraint(m, P_hat[i,t]==r[i]*l[i,t]+P_[i,t]);
        @constraint(m,Q_hat[i,t]==x[i]*l[i,t]+Q_[i,t]);
    end
end

# # #
for i=1:BN
    for t=1:T
        @constraint(m, v[i,t]==v_0[1,t]+(r[i]^2+x[i]^2)*l[i,t]-
            2*(r[i]*P_hat[i,t]+x[i]*Q_hat[i,t]));
    end
end

# #
for i=1:BN
    for t=1:T
        # Power_=v_0[1,t]*l[i,t];
        # @constraint(m, (P_hat[i,t]^2 + Q_hat[i,t]^2 + v_^2 )<= u_^2)
        # @constraint(m, sum([P_hat[i,t]; Q_hat[i,t]; v_].^2 )<= u_^2)
        A_=[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1];
        c_=[0; 0; 1; 0]
        @constraint(m, norm(A_*[P_hat[i,t]; Q_hat[i,t]; u_[i,t]; v_[i,t]])<=
            dot(c_,[P_hat[i,t]; Q_hat[i,t]; u_[i,t]; v_[i,t]]));
        # @constraint(m, (P_hat[i,t]^2+Q_hat[i,t]^2) <= v_0[1,t]^2)
    end
end

for i=1:BN
    for t=1:T
        @constraint(m, u_[i,t]==0.5*(v_0[1,t]+l[i,t]))
        @constraint(m, v_[i,t]==0.5*(v_0[1,t]-l[i,t]))
    end
end
# println(m)
# ################################################################################

@objective(m, Min, fn_cost(delta_t,lambda,P_hat,beta,Pg,Pg_max,alpha,P_rsrv))
status=solve(m);
# println("status =", status)


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

optimized_value=Optimization_output_struct(
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
return optimized_value

end


function fn_cost(delta_t,lambda,P_hat,beta,Pg,Pg_max,alpha,P_rsrv)
    PP3=P_rsrv;
    PP=sum(P_hat,dims=1);
    PP2=sum((Pg-Pg_max));
    return delta_t*(lambda[:,1]'*PP'-beta*PP2-alpha[:,:]'*PP3[:])[1,1]
end

function GML_output(opt_value)
###############################################################################
# write the solar file
println("GML - Start writing files!");
total_n=size(opt_value.Pg)[2];
time_stamp=[1:total_n]
# Vector(Int64, size(opt_value.Pg)[2])
time_stamp=zeros(1,total_n)
for i=1:size(opt_value.Pg)[2]
    time_stamp[i]=i
end

Solar=[time_stamp; opt_value.Pg; opt_value.Qg]';
# println(size(Solar))
CSV.write("solar.csv", DataFrame(Solar, [:timestamp, :Pg_1, :Pg_2, :Pg_3,
    :Pg_4, :Pg_5, :Pg_6, :Pg_7, :Pg_8, :Pg_9, :Pg_10, :Pg_11, :Pg_12, :Qg_1,
    :Qg_2, :Qg_3, :Qg_4, :Qg_5, :Qg_6, :Qg_7, :Qg_8, :Qg_9, :Qg_10,
    :Qg_11, :Qg_12]));

# write the battery file
Battery=[time_stamp; opt_value.R; opt_value.C; opt_value.B]';
# println(size(Battery))
CSV.write("battery.csv", DataFrame(Battery, [:timestamp, :R_1, :R_2, :R_3,
    :R_4, :R_5, :R_6, :R_7, :R_8, :R_9, :R_10, :R_11, :R_12, :C_1,
    :C_2, :C_3, :C_4, :C_5, :C_6, :C_7, :C_8, :C_9, :C_10,
    :C_11, :C_12, :B_1, :B_2, :B_3, :B_4, :B_5, :B_6, :B_7, :B_8, :B_9, :B_10,
     :B_11, :B_12]));

println("GML - Finish writting files! ")
end
