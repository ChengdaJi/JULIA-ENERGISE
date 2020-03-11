function price_traj_det(t, ancillary_type, price_raw, T)
# this function generates the price trajectory
    PTID=3;
    lambda_rt=price_raw.LMP_RT[t+1]
    lambda_scenario = ones(6,1)*reshape(price_raw.LMP_RT[t+2:t+T], 1, T-1);
    if ancillary_type == "without"
        alpha_rt=0;
        alpha_scenario=zeros(6,T-1);
    elseif ancillary_type == "10min"
        alpha_rt=price_raw.RSRV_10[t+1];
        alpha_scenario=ones(6,1)*reshape(price_raw.RSRV_10[t+2:t+T], 1, T-1);
    elseif ancillary_type == "30min"
        alpha_rt=price_raw.RSRV_30[t];
        alpha_scenario=ones(6,1)*reshape(price_raw.RSRV_30[t+2:t+T], 1, T-1);
    end
    probability=[1 0 0 0 0 0];
    price = price_struct(lambda_rt, lambda_scenario, alpha_rt, alpha_scenario, probability)
    return price;
end

function pd_traj_det(t, pd_raw, T)
    # this function generates the demand trajectory
    traj=pd_raw.pd_rt[:, t+1:t+T];
    da=pd_raw.pd_da[:,t+1:t+T];
    sigma = zeros(12, T);
    # pd = pd_struct(traj, sigma);
    pd = (da=(da), traj=(traj), sigma=(sigma));
    return pd
end

function pg_traj_det(t, pg_raw, p_rate, T);
    rt_raw = p_rate*pg_raw.pg_rt;
    da_raw = p_rate*pg_raw.pg_da;
    mu = rt_raw[:,t+1:t+T];
    da = da_raw[:,t+1:t+T];
    mu_ct = rt_raw[:,t+1];
    mu_scenario=rt_raw[:,t+1:t+T-1];
    sigma = zeros(12,T-1)
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=(da=(da),mu=(mu),mu_rt=(mu_rt),mu_scenario=(mu_scenario),sigma=(sigma));
    return pg
end
