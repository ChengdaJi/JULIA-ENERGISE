function GML_optmization(obj)
    println("GML - Start performing optimization");
    #################
    #  F,T or BN, T #
    #################
    ## The parameters

    T = obj.T; # number of time slots
    F = obj.F; # number of feeders
    BN = obj.BN; # number of banks
    Pd = obj.Pd; # demand data (input)
    Qd = obj.Qd; # demand (input)

    Pg_max = obj.Pg_max; # maximum power generation
    Pg_min = obj.Pg_min; # minimum power generation

    Qg_max = obj.Qg_max; # maximum reactive power generation
    Qg_min = obj.Qg_min; # minimum eactive power generation

    B_max = obj.B_max; # maximum storage level
    B_min = obj.B_min; # minimum storage level

    R_max = obj.R_max; # the maximum discharge rate
    R_min = obj.R_min; # the minimum discharge rate

    C_max = obj.C_max; # maximum storage level
    C_min = obj.C_min; # minimum storage level

    b_0 = obj.b_0; # storage level
    W = obj.W; # exogenour change
    delta_t = obj.delta_t; # time interval

    P_rsrv_min = obj.P_rsrv_min; # The minimum ancillary power

    S=obj.S; # the apparents power

    V_min=obj.V_min;
    V_max=obj.V_max;
    alpha=obj.alpha; # price in ancillary power
    beta=obj.beta; # price in cost function
    lambda= obj.lambda; # the grid spot price
    tau = obj.tau;

    r=obj.r; # the resistance
    x=obj.x; # the reactance
    k=obj.k; # time with ancellary

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

optimied_value=Optmization_output(
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
return optimied_value

end


function fn_cost(delta_t,lambda,P_hat,beta,Pg,Pg_max,alpha,P_rsrv)
    PP3=P_rsrv;
    PP=sum(P_hat,dims=1);
    PP2=sum((Pg-Pg_max));
    return delta_t*(lambda[:,1]'*PP'+beta*PP2-alpha[:,:]'*PP3[:])[1,1]
end
