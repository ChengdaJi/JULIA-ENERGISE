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
    price = price_struct(lambda_rt, lambda_scenario, alpha_rt, alpha_scenario, probability)
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
    # sigma_1 = hcat(pred.*(ones(12,1)*reshape(pd_MAPE["demand_MAPE_2hr"][1,1:Pred_length].^2, 1, Pred_length)),
    #     zeros(12, T-1-Pred_length));
    ratio = reshape(0.01:0.01/23:0.01/23*(Pred_length-1)+0.01, 1, Pred_length)
    sigma_1 = hcat(ones(12,1)*ratio.*pred, zeros(12, T-1-Pred_length))
    traj_1 = hcat(rt, pred, pd_raw.pd_da[:,t+1+Pred_length:t+T-1]);
    sigma_c=[];
    traj_c=[]
    for i=1:BN/3
        push!(sigma_c,sigma_1);
        push!(traj_c,traj_1);
    end
    traj=vcat(traj_c...)
    sigma=vcat(sigma_c...)
    # pd = pd_struct(traj, sigma);
    pd = pd_struct(traj, sigma);
    return pd
end

function pg_traj(t, pg_raw, pg_noise, solar_error_max, p_rate, T, Pred_length);
    rt_raw = p_rate*pg_raw.pg_rt;
    da_raw = p_rate*pg_raw.pg_da;
    mu_rt = rt_raw[:,t]
    pred = zeros(12, Pred_length)
    for feeder=1:12
        pg_ratio = positive_array(pg_noise[feeder][t,1:Pred_length].+1);
        pred[feeder, :]=rt_raw[feeder,t+1:t+Pred_length].*pg_ratio;
    end
    mu_scenario = hcat(pred, da_raw[:,t+1+Pred_length:t+T-1])
    mu = hcat(reshape(mu_rt, 12, 1), mu_scenario);
    ratio = reshape(0.01:(solar_error_max-0.01)/23:(solar_error_max-0.01)/23*(Pred_length-1)+0.01, 1, Pred_length)
    sigma = hcat(ones(12,1)*ratio.*pred, zeros(12, T-1-Pred_length));
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    return pg
end
