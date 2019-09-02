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
    probability=delta_rt_raw["probability"][PTID,t];
    price = price_struct(lambda_rt, lambda_scenario, 0, 0, probability)
    return price;
end

# function pd_traj(t, pd_raw_one, pd_MAPE, e, BN,T, Pred_length)
#     # this function generates the demand trajectory
#     pd_raw = hcat(pd_raw_one["Pd"], 0.9*pd_raw_one["Pd"])/1000;
#     rt=pd_raw[:,t];
#     pred=pd_raw[:,t+1:t+Pred_length];
#     sigma_1 = hcat(pred.*(ones(12,1)*reshape(pd_MAPE["demand_MAPE_2hr"][1,1:Pred_length].^2, 1, Pred_length)),
#         zeros(12, T-1-Pred_length));
#     traj_1 = hcat(rt, pd_raw[:,t+1:t+T-1]);
#     sigma_c=[];
#     traj_c=[]
#     for i=1:BN/3
#         push!(sigma_c,sigma_1);
#         push!(traj_c,traj_1);
#     end
#     traj=vcat(traj_c...)
#     sigma=vcat(sigma_c...)
#     # pd = pd_struct(traj, sigma);
#     pd = pd_struct(traj, sigma);
#     return pd
# end

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


# function pd_traj(t, pd_raw, pd_MAPE, e, BN,T, Pred_length)
#     # this function generates the demand trajectory
#     rt=pd_raw[:,t];
#     pred=pd_raw[:,t+1:t+Pred_length];
#     # sigma_1 = hcat(pred.*(ones(12,1)*reshape(pd_MAPE["demand_MAPE_2hr"][1,1:Pred_length].^2, 1, Pred_length)),
#     #     zeros(12, T-1-Pred_length));
#     sigma_1 = zeros(12, 287)
#     traj_1 = hcat(rt, pd_raw[:,t+1:t+T-1]);
#     sigma_c=[];
#     traj_c=[]
#     for i=1:BN/3
#         push!(sigma_c,sigma_1);
#         push!(traj_c,traj_1);
#     end
#     traj=vcat(traj_c...)
#     sigma=vcat(sigma_c...)
#     # pd = pd_struct(traj, sigma);
#     pd = pd_struct(traj, sigma);
#     return pd
# end

# function pg_traj(t, pg_raw_one, pg_RMSE, error, pd_traj, p_rate, BN,T, Pred_length)
#     pg_traj1_raw=hcat(pg_raw_one["Pg1"], pg_raw_one["Pg1"]);
#     pd_current_sum_1=sum(pd_traj[1:8,1:T]);
#     pd_portion1=zeros(8,1)
#     for feeder=1:8
#         pd_portion1[feeder] = sum(pd_traj[feeder,:])/pd_current_sum_1;
#     end
#     mu1=pd_portion1*reshape(pg_traj1_raw[t:t+T-1],1,T)/sum(pg_raw_one["Pg1"])*pd_current_sum_1*p_rate;
#
#     pg_traj2_raw=hcat(pg_raw_one["Pg2"], pg_raw_one["Pg2"]);
#     pd_current_sum_2=sum(pd_traj[9:12,:]);
#     pd_portion2=zeros(4,1)
#     for feeder=1:4
#         pd_portion2[feeder] = sum(pd_traj[feeder,:])/pd_current_sum_2;
#     end
#     mu2=pd_portion2*reshape(pg_traj2_raw[t:t+T-1],1,T)/sum(pg_raw_one["Pg2"])*pd_current_sum_2*p_rate;
#     mu_1=vcat(mu1,mu2);
#     mu_rt_1=reshape(mu_1[:,1],12,1)
#     mu_scenario_1=zeros(12,T-1);
#     sigma_1 = zeros(12,T-1)
#     # for feeder=1:12
#     #     for time=2:T
#     #         if mu_1[feeder, time]==0
#     #         else
#     #             mu_scenario_1[feeder, time-1] = positive_array(1+error[current_time, time]*mu_1[feeder, time])
#     #             sigma_1[feeder, time-1] =
#     #     end
#     # end
#
#
#     sigma_1=hcat((mu_1[:,2:Pred_length+1].*(ones(12,1)*reshape(pg_RMSE["solar_RRMSE_2hr"][1,1:Pred_length],1,Pred_length))).^2,zeros(12,T-1-Pred_length));
#     sigma_c=[];
#     mu_c=[];
#     mu_rt_c=[];
#     mu_scenario_c=[];
#     for i=1:BN/3
#         push!(sigma_c,sigma_1);
#         push!(mu_c,mu_1);
#         push!(mu_rt_c,mu_rt_1);
#         push!(mu_scenario_c,mu_scenario_1);
#     end
#     mu=vcat(mu_c...)
#     mu_rt=vcat(mu_rt_c...)
#     mu_scenario=vcat(mu_scenario_c...)
#     sigma=vcat(sigma_c...)
#     # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
#     pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
#     return pg
# end

# function pg_traj(t, pg_raw_one, pg_noise, pd_traj, p_rate, BN,T, Pred_length);
#     pg_traj1_raw=hcat(pg_raw_one["Pg1"], pg_raw_one["Pg1"]);
#     pd_current_sum_1=sum(pd_traj[1:8,1:T]);
#     pd_portion1=zeros(8,1)
#     for feeder=1:8
#         pd_portion1[feeder] = sum(pd_traj[feeder,1:T])/pd_current_sum_1;
#     end
#     println(sum(pd_portion1))
#     mu1=pd_portion1*reshape(pg_traj1_raw[t:t+T-1],1,T)/sum(pg_raw_one["Pg1"])*pd_current_sum_1*p_rate;
#
#     pg_traj2_raw=hcat(pg_raw_one["Pg2"], pg_raw_one["Pg2"]);
#     pd_current_sum_2=sum(pd_traj[9:12,:]);
#     pd_portion2=zeros(4,1)
#     for feeder=1:4
#         pd_portion2[feeder] = sum(pd_traj[feeder+8,1:288])/pd_current_sum_2;
#     end
#     println(sum(pd_portion2))
#     mu2=pd_portion2*reshape(pg_traj2_raw[t:t+T-1],1,T)/sum(pg_raw_one["Pg2"])*pd_current_sum_2*p_rate;
#     mu_1=vcat(mu1,mu2);
#     mu_rt_1=reshape(mu_1[:,1],12,1)
#     mu_scenario_1=zeros(12,T-1);
#     sigma_1 = zeros(12,T-1)
#     # for feeder=1:12
#     #     for time=2:T
#     #         if mu_1[feeder, time]==0
#     #         else
#     #             mu_scenario_1[feeder, time-1] = positive_array(1+error[current_time, time]*mu_1[feeder, time])
#     #             sigma_1[feeder, time-1] =
#     #     end
#     # end
#
#
#     # sigma_1=hcat((mu_1[:,2:Pred_length+1].*(ones(12,1)*reshape(pg_RMSE["solar_RRMSE_2hr"][1,1:Pred_length],1,Pred_length))).^2,zeros(12,T-1-Pred_length));
#     sigma_1 = zeros(12,287)
#     sigma_c=[];
#     mu_c=[];
#     mu_rt_c=[];
#     mu_scenario_c=[];
#     for i=1:BN/3
#         push!(sigma_c,sigma_1);
#         push!(mu_c,mu_1);
#         push!(mu_rt_c,mu_rt_1);
#         push!(mu_scenario_c,mu_scenario_1);
#     end
#     mu=vcat(mu_c...)
#     mu_rt=vcat(mu_rt_c...)
#     mu_scenario=vcat(mu_scenario_c...)
#     sigma=vcat(sigma_c...)
#     # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
#     pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
#     return pg
# end



function pg_traj(t, pg_raw, pg_noise, solar_error_max, p_rate, T, Pred_length);
    rt_raw = p_rate/0.5*pg_raw.pg_rt;
    da_raw = p_rate/0.5*pg_raw.pg_da;
    println(size(rt_raw))
    println(size(da_raw))
    mu_rt = rt_raw[:,t]
    println(size(mu_rt))
    pred = zeros(12, Pred_length)
    for feeder=1:12
        pg_ratio = positive_array(pg_noise[feeder][t,1:Pred_length].+1);
        pred[feeder, :]=rt_raw[feeder,t+1:t+Pred_length].*pg_ratio;
    end
    println(size(pred))
    println(size(da_raw[:,t+1+Pred_length:t+T-1]))
    mu_scenario = hcat(pred, da_raw[:,t+1+Pred_length:t+T-1])
    mu = hcat(reshape(mu_rt, 12, 1), mu_scenario);
    println(size(mu))
    ratio = reshape(0.01:(solar_error_max-0.01)/23:(solar_error_max-0.01)/23*(Pred_length-1)+0.01, 1, Pred_length)
    println(size(ratio))
    sigma = hcat(ones(12,1)*ratio.*pred, zeros(12, T-1-Pred_length));
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    return pg
end
