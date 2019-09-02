struct pd_struct
    # trajectory
    traj;
    # variance
    sigma;
end

struct price_struct
    # energy market price
        # real_time price
        lambda_rt;
        # offline price
        lambda_scenario;
    # ancillary price
        # real_time ancillary
        alpha_rt;
        # offline price
        alpha_scenario;
    # probability
        probability;
end

struct pg_struct
    # average
    mu;
    # average realtime
    mu_rt;
    # average scenario
    mu_scenario;
    # sigma
    sigma;
end

struct bd_struct
        #number of time slots
        T;
        #number of feeders
        F;
        #number of banks
        BN;
        # feeder level (Solar virtual storage and reactive power)
        # minimun power output
        Pg_min;
        # maximum reactive power output
        Qf_max;
        # minimun reactive power output
        Qf_min;
        # the manimum power of storage
        B_max;
        # the minimum power of storage
        B_min;
        # the maximum power discharge rate
        R_max;
        # the minimum power discharge rate
        R_min;
        # the exogenour change
        W;
        # time interval
        delta_t;
        ## Ancillary Services
        # The minimum ancillary power
        P_rsrv_min;
        # reserve power delay
        tau;
        # reserve power serving time
        k;
        # radial graph constraints
        # the maximum voltage
        V_max;
        # the minimum voltage
        V_min;
        # impedance
        r;
        # impedance
        x;
        # penalty on solar
        beta;
        # apperatne power
        S;
    end

struct Optimization_output_struct
    # power at feeder
    Pf;
    # reactive power at feeder
    Qf;
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
