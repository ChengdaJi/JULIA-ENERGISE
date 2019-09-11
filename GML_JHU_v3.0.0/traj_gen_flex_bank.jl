function price_traj(t, ancillary_type, price_raw, delta_rt_raw, T, Pred_length)
# this function generates the price trajectory
    PTID=3;
    lambda_rt=price_raw.LMP_RT[t]
    # lambda_rt=price_raw["lambda_rt"][t]
    lambda_pd=lambda_rt.+delta_rt_raw["RTcentroid"][PTID,t][:,1:Pred_length];
    Pred_scen=length(lambda_pd[:,1]);

    # println(size(lambda_pd))
    lambda_da=price_raw.LMP_DA[t+Pred_length+1:t+T-1];
    # lambda_da=price_raw["lambda_da"][t+Pred_length+1:t+T-1];
    # println(size(lambda_da))
    lambda_offline=ones(Pred_scen,1)*reshape(lambda_da, 1, length(lambda_da));
    lambda_scenario = hcat(lambda_pd,lambda_offline);
    lambda_hr_oneset = zeros(1, 24-ceil(Int,T/12))
    for hr=1:21
        lambda_hr_oneset[hr] = price_raw.LMP_DA[t+T+12*(hr-1)]
    end
    lambda_hr = ones(Pred_scen,1)*lambda_hr_oneset;
    if ancillary_type == "without"
        alpha_rt=0;
        alpha_scenario=zeros(Pred_scen,T-1);
    elseif ancillary_type == "10min"
        alpha_rt=price_raw.RSRV_10[t];
        alpha_pd=alpha_rt.+delta_rt_raw["RSRV10centroid"][PTID,t][:,1:Pred_length];
        alpha_scenario=hcat(alpha_pd,zeros(Pred_scen,T-1-Pred_length));
    elseif ancillary_type == "30min"
        alpha_rt=price_raw.RSRV_30[t];
        alpha_pd=alpha_rt.+delta_rt_raw["RSRV30centroid"][PTID,t][:,1:Pred_length];
        alpha_scenario=hcat(alpha_pd,zeros(Pred_scen,T-1-Pred_length));
    end
    probability=delta_rt_raw["probability"][PTID,t];
    price = (lambda_rt=(lambda_rt), lambda_scenario=(lambda_scenario),alpha_rt=(alpha_rt),
        alpha_scenario=(alpha_scenario), probability=(probability), lambda_hr=(lambda_hr))
    return price;
end

function pd_traj(t, pd_raw, pd_noise, BN,T, Pred_length)
    # this function generates the demand trajectory
    rt=pd_raw.pd_rt[:,t];
    pred = zeros(12, Pred_length)
    for feeder=1:12
        pd_ratio = positive_array(pd_noise[feeder][t,1:Pred_length].+1);
        pred[feeder, :]=pd_raw.pd_rt[feeder,t+1:t+Pred_length].*pd_ratio;
    end
    pred_hr_oneset = zeros(12,24-ceil(Int,T/12));
    for hr=1:21
        pred_hr_oneset[:,hr] = pd_raw.pd_da[:,t+T+12*(hr-1)]
    end
    ratio = reshape(0.01:0.01/23:0.01/23*(Pred_length-1)+0.01, 1, Pred_length)
    sigma = repeat(hcat(ones(12,1)*ratio.*pred, zeros(12, T-1-Pred_length)), trunc(Int, BN/3),1);
    traj = repeat(hcat(rt, pred, pd_raw.pd_da[:,t+1+Pred_length:t+T-1]),trunc(Int, BN/3),1);
    pred_hr = repeat(pred_hr_oneset, trunc(Int, BN/3),1);
    pd = (traj=(traj), sigma=(sigma), pred_hr=(pred_hr))
    return pd
end

function pg_traj(t, pg_raw, pg_noise, solar_error_max, p_rate, BN, T, Pred_length);
    rt_raw = p_rate*pg_raw.pg_rt;
    da_raw = p_rate*pg_raw.pg_da;
    mu_rt = repeat(reshape(rt_raw[:,t],12,1),trunc(Int, BN/3),1);
    pred = zeros(12, Pred_length)
    for feeder=1:12
        pg_ratio = positive_array(pg_noise[feeder][t,1:Pred_length].+1);
        pred[feeder, :]=rt_raw[feeder,t+1:t+Pred_length].*pg_ratio;
    end
    mu_scenario = repeat(hcat(pred, da_raw[:,t+1+Pred_length:t+T-1]),trunc(Int, BN/3),1);
    mu = hcat(mu_rt, mu_scenario);
    ratio = reshape(0.01:(solar_error_max-0.01)/23:(solar_error_max-0.01)/23*(Pred_length-1)+0.01, 1, Pred_length)
    sigma = repeat(hcat(ones(12,1)*ratio.*pred, zeros(12, T-1-Pred_length)),trunc(Int, BN/3),1);

    pred_hr_oneset = zeros(12,24-ceil(Int,T/12))
    for hr=1:21
        pred_hr_oneset[:,hr] = da_raw[:,t+T+12*(hr-1)]
    end
    pred_hr = repeat(pred_hr_oneset, trunc(Int, BN/3),1);
    pg=(mu=(mu),mu_rt=(mu_rt),mu_scenario=(mu_scenario),sigma=(sigma), pred_hr=(pred_hr));
    return pg
end
