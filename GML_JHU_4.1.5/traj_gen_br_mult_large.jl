function price_traj(t, ancillary_type, price_raw, delta_rt_raw, T, Pred_length)
# this function generates the price trajectory
    PTID=3;
    lambda_rt=price_raw.LMP_RT[t]
    lambda_ct=price_raw.LMP_RT[t+1]
    # lambda_rt=price_raw["lambda_rt"][t]
    lambda_pd=lambda_rt.+delta_rt_raw["RTcentroid"][PTID,t][:,1:Pred_length];
    Pred_scen=length(lambda_pd[:,1]);

    # println(size(lambda_pd))
    lambda_da=price_raw.LMP_DA[t+Pred_length+1:t+T];
    # lambda_da=price_raw["lambda_da"][t+Pred_length+1:t+T-1];
    # println(size(lambda_da))
    lambda_offline=ones(Pred_scen,1)*reshape(lambda_da, 1, length(lambda_da));
    lambda_scenario = hcat(lambda_pd,lambda_offline);
    if ancillary_type == "without"
        alpha_ct=0;
        alpha_scenario=zeros(Pred_scen,T);
    elseif ancillary_type == "10min"
        alpha_rt=price_raw.RSRV_10[t];
        alpha_ct=price_raw.RSRV_10[t+1];
        # for t ==  4
        # elseif t == 
        alpha_pd=alpha_rt.+delta_rt_raw["RSRV10centroid"][PTID,t][:,1:Pred_length];
        alpha_scenario=hcat(alpha_pd, zeros(Pred_scen,T-Pred_length));
        # alpha_scenario = 50*ones(Pred_scen, T)
        # println(sum(alpha_scenario))
    elseif ancillary_type == "30min"
        alpha_rt=price_raw.RSRV_30[t];
        alpha_ct=price_raw.RSRV_30[t+1];
        alpha_pd=alpha_rt.+delta_rt_raw["RSRV30centroid"][PTID,t][:,1:Pred_length];
        alpha_scenario=hcat(alpha_pd,zeros(Pred_scen,T-Pred_length));
    end
    probability=delta_rt_raw["probability"][PTID,t];

    price = (lambda_ct=(lambda_ct), lambda_scenario = (lambda_scenario),
        alpha_ct=(alpha_ct), alpha_scenario=(alpha_scenario),
        probability=(probability))

    return price;
end

function pd_traj_large(t, pd_raw, pd_noise,T, Pred_length, base, mult)
    # println(mult)
    # this function generates the demand trajectory
    rt=pd_raw.pd_rt[:,t];
    pred = zeros(12*mult, Pred_length)
    for feeder=1:12*mult
        # println(feeder)
        feeder_noise = Int(feeder-12*floor((feeder-1)/12))
        # println(feeder_noise)
        pd_ratio = positive_array(pd_noise[feeder_noise][t,1:Pred_length].+1);
        pred[feeder, :]=pd_raw.pd_rt[feeder,t+1:t+Pred_length].*pd_ratio;
    end

    da = pd_raw.pd_da[:,t+1+Pred_length:t+T];
    ratio = reshape(0.01:0.01/23:0.01/23*(Pred_length-1)+0.01, 1, Pred_length)
    pred_square = pred.^2;

    traj_bus = zeros(15*mult, 288);
    traj_feeder = hcat(pred, da);
    sigma_bus = zeros(15*mult, 288);
    sigma_feeder = hcat(ones(12*mult,1)*ratio, 0.08*da.^2);
    ct_bus = zeros(15*mult, 1);
    ct_feeder = pd_raw.pd_rt[:,t+1];
    for set = 1:mult
        traj_bus[(set-1)*15+1:set*15,:]=vcat(
            zeros(3,288),traj_feeder[(set-1)*12+1:set*12, :])
        sigma_bus[(set-1)*15+1:set*15,:]=vcat(
            zeros(3,288),sigma_feeder[(set-1)*12+1:set*12, :])
        ct_bus[(set-1)*15+1:set*15]=vcat(
            zeros(3,1), ct_feeder[(set-1)*12+1:set*12])
    end

    sigma = vcat(zeros(1,288), sigma_bus)
    traj = vcat(zeros(1,288), traj_bus)
    ct=vcat(0, ct_bus);

    pd = (traj=(traj/base.MVA), sigma = (sigma/(base.MVA)/base.MVA),
    ct=(ct/base.MVA))
    return pd
end

function pg_traj_large(t, pg_raw, pg_noise, solar_error_max, p_rate, T,
    Pred_length, base, mult)
    rt_raw = p_rate*pg_raw.pg_rt;
    da_raw = p_rate*pg_raw.pg_da;
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
