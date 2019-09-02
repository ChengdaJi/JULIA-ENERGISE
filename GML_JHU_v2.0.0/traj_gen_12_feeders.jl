function price_traj(t, ancillary_type, price_raw, delta_rt_raw)
# this function generates the price trajectory
    PTID=3;
    lambda_rt=price_raw["lambda_rt"][t]
    lambda_pd=lambda_rt.+delta_rt_raw["RTcentroid"][PTID,t];
    Pred_length=length(lambda_pd[1,:]);
    Pred_scen=length(lambda_pd[:,1]);
    # println(size(lambda_pd))
    lambda_da=price_raw["lambda_da"][t+Pred_length+1:t+288-1];
    # println(size(lambda_da))
    lambda_offline=ones(Pred_scen,1)*reshape(lambda_da, 1, length(lambda_da));
    lambda_scenario = hcat(lambda_pd,lambda_offline);
    # println(size(lambda_scenario))
    if ancillary_type == "without"
        alpha_rt=0;
        alpha_scenario=zeros(6,287);
    elseif ancillary_type == "10min"
        alpha_rt=price_raw["alpha_rt_10"][t];
        alpha_pd=alpha_rt.+delta_rt_raw["RSRV10centroid"][PTID,t];
        alpha_scenario=hcat(alpha_pd,zeros(6,287-24));
    elseif ancillary_type == "30min"
        alpha_rt=price_raw["alpha_rt_30"][t];
        alpha_pd=alpha_rt.+delta_rt_raw["RSRV30centroid"][PTID,t];
        alpha_scenario=hcat(alpha_pd,zeros(6,287-24));
    end
    probability=delta_rt_raw["probability"][PTID,t];
    price = price_struct(lambda_rt, lambda_scenario, alpha_rt, alpha_scenario, probability)
    return price;
end

function pd_traj(t, pd_raw_one, pd_MAPE, e,BN)
    # this function generates the demand trajectory
    pd_raw = hcat(pd_raw_one["Pd"], 0.9*pd_raw_one["Pd"])/1000;
    rt=pd_raw[:,t]+e;
    pred=pd_raw[:,t+1:t+24];
    # testthing = (ones(12,1)*reshape(pd_MAPE["demand_MAPE_2hr"]),1,24).^2
    sigma_1 = hcat(pred.*(ones(12,1)*pd_MAPE["demand_MAPE_2hr"]).^2, zeros(12, 287-25+1));
    traj_1 = hcat(rt, pd_raw[:,t+1:t+287]);
    sigma_c=[];
    traj_c=[]
    for i=1:BN/3
        push!(sigma_c,sigma_1);
        push!(traj_c,traj_1);
    end
    traj=vcat(traj_c...)
    sigma=vcat(sigma_c...)
    pd = pd_struct(traj, sigma);
    return pd
end

function pg_traj(t, pg_raw_one, pg_RMSE, pd_traj, p_rate,BN)
    pg_traj1_raw=hcat(pg_raw_one["Pg1"], pg_raw_one["Pg1"]);
    pd_current_sum_1=sum(pd_traj[1:8,:]);
    pd_portion1=zeros(8,1)
    for feeder=1:8
        pd_portion1[feeder] = sum(pd_traj[feeder,:])/pd_current_sum_1;
    end
    mu1=pd_portion1*reshape(pg_traj1_raw[t:t+287],1,288)/sum(pg_raw_one["Pg1"])*pd_current_sum_1*p_rate;

    pg_traj2_raw=hcat(pg_raw_one["Pg2"], pg_raw_one["Pg2"]);
    pd_current_sum_2=sum(pd_traj[9:12,:]);
    pd_portion2=zeros(4,1)
    for feeder=1:4
        pd_portion2[feeder] = sum(pd_traj[feeder,:])/pd_current_sum_2;
    end
    mu2=pd_portion2*reshape(pg_traj2_raw[t:t+287],1,288)/sum(pg_raw_one["Pg2"])*pd_current_sum_2*p_rate;
    mu_1=vcat(mu1,mu2);
    mu_rt_1=reshape(mu_1[:,1],12,1)
    mu_scenario_1=mu_1[:,2:end];
    sigma_1=hcat((mu_1[:,2:25].*(ones(12,1)*reshape(pg_RMSE["solar_RRMSE_2hr"],1,24))).^2,zeros(12,287-24));
    sigma_c=[];
    mu_c=[];
    mu_rt_c=[]
    mu_scenario_c=[]
    for i=1:BN/3
        push!(sigma_c,sigma_1);
        push!(mu_c,mu_1);
        push!(mu_rt_c,mu_rt_1);
        push!(mu_scenario_c,mu_scenario_1);
    end
    mu=vcat(mu_c...)
    mu_rt=vcat(mu_rt_c...)
    mu_scenario=vcat(mu_scenario_c...)
    sigma=vcat(sigma_c...)
    pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    return pg
end
