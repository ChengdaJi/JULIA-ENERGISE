struct demand_struct
    Pd
    Qd
end

struct price_struct
    alpha
    beta
    lambda
end

struct Pg_trend_struct
    Pg1
    Pg2
end

struct bd_struct
        #number of time slots
        T;
        #number of feeders
        F;
        #number of banks
        BN;
        # demand
        Pd;
        # reactive demand
        Qd;
        # maximum power output
        Pg_max;
        # minimun power output
        Pg_min;
        # maximum reactive power output
        Qg_max;
        # minimun reactive power output
        Qg_min;
        # the manimum power of storage
        B_max;
        # the minimum power of storage
        B_min;
        # the maximum power discharge rate
        R_max;
        # the minimum power discharge rate
        R_min;
        # power storage
        C_max;
        # the minimum power of storage
        C_min;
        # power storage
        b_0;
        # the exogenour change
        W;
        # time interval
        delta_t;
        # The maximum ancillary power
        P_rsrv_min;
        # the apparents power
        S;
        # the maximum voltage
        V_max;
        # the minimum voltage
        V_min;
        # price of ancillary power
        alpha;
        # price of cost
        beta;
        # the grid spot price
        lambda;
        tau;
        # impedance
        r;
        # impedance
        x;
        # time asso/w ancellary
        k;
    end

struct Optimization_output_struct
    # solar generation
    Pg;
    # reactive solar generation
    Qg;
    # battery SOC
    B;
    # battery power discharge
    R;
    # battery reactive power discharge
    C;
    # branch power flow
    P_hat;
    # branch reactive power flow
    Q_hat;
    # current square
    l;
    # voltage square
    v;
    # schedule ancillary power
    P_rsrv;
    # schedule ancillary energy
    B_rsrv;
    # power consuming at head node
    P_0;
    # operating cost
    cost;
end
