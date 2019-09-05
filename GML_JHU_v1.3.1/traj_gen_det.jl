function price_traj_det(t, ancillary_type, price_raw, T)
# this function generates the price trajectory
    PTID=3;
    lambda_rt=price_raw.LMP_RT[t]
    lambda_scenario = ones(6,1)*reshape(price_raw.LMP_RT[t+1:t+T-1], 1, T-1);
    # println(size(lambda_scenario))
    # if ancillary_type == "without"
    #     alpha_rt=0;
    #     alpha_scenario=zeros(Pred_scen,T-1);
    # elseif ancillary_type == "10min"
    #     alpha_rt=price_raw["alpha_rt_10"][t];
    #     alpha_pd=alpha_rt.+delta_rt_raw["RSRV10centroid"][PTID,t][:,1:Pred_length];
    #     alpha_scenario=hcat(alpha_pd,zeros(Pred_scen,T-1-Pred_length));
    # elseif ancillary_type == "30min"
    #     alpha_rt=price_raw["alpha_rt_30"][t];
    #     alpha_pd=alpha_rt.+delta_rt_raw["RSRV30centroid"][PTID,t][:,1:Pred_length];
    #     alpha_scenario=hcat(alpha_pd,zeros(Pred_scen,T-1-Pred_length));
    # end
    probability=[1, 0, 0, 0, 0, 0];
    price = price_struct(lambda_rt, lambda_scenario, 0, 0, probability)
    return price;
end

function pd_traj(t, pd_raw, T)
    # this function generates the demand trajectory
    traj=pd_raw.pd_rt[:, t:t+T-1];
    sigma = zeros(12, T-1);
    # pd = pd_struct(traj, sigma);
    pd = pd_struct(traj, sigma);
    return pd
end

function pg_traj(t, pg_raw, T);
    rt_raw = pg_raw.pg_rt/23*8.61;
    mu = rt_raw[:,t:t+T-1];
    mu_rt = rt_raw[:,t];
    mu_scenario=rt_raw[:,t+1:t+T-1];
    sigma = zeros(12,T-1)
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    return pg
end
