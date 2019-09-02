function GML_boundaries(Pd, Pg_trend, Lambda, Beta, Alpha)
    println("GML - Start defining boundaries");
###############################################################################
    # demand
    # Pd=demand.Pd/1000
    number_=size(Pd)
    F=number_[1];
    T=number_[2];

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
    Pg_trends1=Pg_trends[1,:]';
    Pg_trends2=Pg_trends[2,:]';
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



    obj_return=obj_struct(T, F, 3, Pd, Qd, Pg_max, Pg_min, Qg_max, Qg_min, B_max,
            B_min, R_max, R_min, C_max, C_min, b_0, W, delta_t,
            P_rsrv_min,S,V_max, V_min, Alpha, Beta, Lambda, tau, r, x, k);
    println("GML - Finish defining boundaries")
    return obj_return

end
