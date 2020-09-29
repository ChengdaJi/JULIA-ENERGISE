function GML_Sys_Ava_large_emergency(T, BN, SN, pd, ancillary_type, icdf,
    B_cap, base, network,pcent_loss)


    mult = Int((BN-1)/15);
    println("===== GML - Boundaries Buildup");
    ###############################################################################
    # feeder level
    Qf_max=0.05*positive_array(pd.traj);
    Qf_min = -Qf_max;
    # minimum solar
    Pg_min = zeros(BN, T);
    # battery
    B_rate=B_cap;
    R_rate=1/3;
    B_max = vcat(zeros(1,T),
        repeat(vcat(zeros(3,T), ones(12,1)*B_rate*ones(1,T)/12), mult));
    println(sum(B_max))
    B_max[5,:] = B_rate*ones(1,T)/12*(1-pcent_loss)

    B_min = zeros(BN,T);
    R_max = R_rate*B_max/base.MVA;
    R_min = -R_max;

    W=zeros(BN,T);
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
    S=35/base.MVA;
    # V
    # Base_V=69;
    # V_min = 0.96;
    # V_max = 1.06;
    bus = network.bus
    branch = network.branch
    V_max = bus[:,12];
    # println(V_max)
    V_min = bus[:,13];
    # impedance
    r=branch[:,3];
    x=branch[:,4];
    # println(r)

    # penalty
    beta=-15;
    println("===== GML - finish building boundaries =====")
    obj=obj_struct(T,0,BN,SN, Pg_min,Qf_max,Qf_min,B_max,B_min,R_max,R_min,W,
    delta_t,P_rsrv_min,tau,k,V_max,V_min,r,x,beta,S,icdf);
    return obj
end

function pg_traj_large_emergency(t, pg_raw, pg_noise, solar_error_max, p_rate, T,
    Pred_length, base, mult, pcent_loss)
    rt_raw = p_rate*pg_raw.pg_rt;
    rt_raw[5,:] = p_rate*pg_raw.pg_rt[5,:]*(1-pcent_loss);
    da_raw = p_rate*pg_raw.pg_da;
    da_raw[5,:] = p_rate*pg_raw.pg_da[5,:]*(1-pcent_loss);

    mu_rt = rt_raw[:,t];
    mu_ct = vcat(zeros(4,1),rt_raw[:,t+1]);
    pred = zeros(12*mult, Pred_length)
    traj025 = 0.01:(0.025-0.01)/(Pred_length-1):0.025
    traj = traj025.+(solar_error_max-0.01)
    for feeder=1:12*mult

        feeder_noise = Int(feeder-12*floor((feeder-1)/12))
        temp=sqrt.(traj./traj025);
        temp1=temp.*pg_noise[feeder_noise][t,1:Pred_length];
        # pg_ratio = positive_array(pg_noise[feeder][t,1:Pred_length].+1);
        pg_ratio=positive_array(temp1.+1)
        pred[feeder, :]=rt_raw[feeder,t+1:t+Pred_length].*pg_ratio;
    end

    da = da_raw[:,t+1+Pred_length:t+T];
    ratio = reshape(traj025.+(solar_error_max-0.01),1,Pred_length)
    pred_square = pred.^2;

    mu_bus = zeros(15*mult, 288);
    mu_feeder = hcat(pred, da);
    sigma_bus = zeros(15*mult, 288);
    sigma_feeder = hcat(
        (ones(12*mult,1)*ratio).*pred_square, p_rate^2*2*solar_error_max*da.^2);
    # sigma_feeder = hcat(
    #     ones(12*mult,1)*ratio.*pred_square,
    #     p_rate^2*2*solar_error_max*ones(12*mult,T-Pred_length));
    mu_ct_bus = zeros(15*mult, 1);
    mu_ct_feeder = rt_raw[:,t+1];
    for set = 1:mult
        mu_bus[(set-1)*15+1:set*15,:]=vcat(
            zeros(3,288),mu_feeder[(set-1)*12+1:set*12, :])
        sigma_bus[(set-1)*15+1:set*15,:]=vcat(
            zeros(3,288),sigma_feeder[(set-1)*12+1:set*12, :])
        mu_ct_bus[(set-1)*15+1:set*15]=vcat(
            zeros(3,1), mu_ct_feeder[(set-1)*12+1:set*12])
    end
    sigma = vcat(zeros(1,288), sigma_bus)
    mu = vcat(zeros(1,288), mu_bus)
    mu_ct=vcat(0, mu_ct_bus);
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=(mu=(mu/base.MVA),sigma=(sigma/(base.MVA)^2), mu_ct=(mu_ct/base.MVA));
    return pg
end
