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
        alpha_pd=alpha_rt.+delta_rt_raw["RSRV10centroid"][PTID,t][:,1:Pred_length];
        alpha_scenario=hcat(alpha_pd,zeros(Pred_scen,T-Pred_length));
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
    # price_struct(lambda_rt, lambda_scenario, alpha_rt, alpha_scenario, probability)
    return price;
end

function price_traj_det(t, ancillary_type, price_raw, T)
# this function generates the price trajectory
    PTID=3;
    lambda_ct=price_raw.LMP_RT[t+1]
    lambda_scenario = ones(6,1)*reshape(price_raw.LMP_RT[t+1:t+T], 1, T);
    if ancillary_type == "without"
        alpha_ct=0;
        alpha_scenario=zeros(6,T);
    elseif ancillary_type == "10min"
        alpha_ct=price_raw.RSRV_10[t+1];
        # println(alpha_ct)
        alpha_scenario=ones(6,1)*reshape(price_raw.RSRV_10[t+1:t+T], 1, T);
    elseif ancillary_type == "30min"
        alpha_ct=price_raw.RSRV_30[t+1];
        alpha_scenario=ones(6,1)*reshape(price_raw.RSRV_30[t+1:t+T], 1, T);
    end
    probability=[1 0 0 0 0 0];
    price = (lambda_ct=(lambda_ct), lambda_scenario = (lambda_scenario),
        alpha_ct=(alpha_ct), alpha_scenario=(alpha_scenario),
        probability=(probability))
    return price;
end

function pd_traj(t, pd_raw, pd_noise, BN,T, Pred_length)
    # this function generates the demand trajectory
    rt=pd_raw.pd_rt[:,t];
    ct=pd_raw.pd_rt[:,t+1];
    pred = zeros(12, Pred_length)
    for feeder=1:12
        pd_ratio = positive_array(pd_noise[feeder][t,1:Pred_length].+1);
        pred[feeder, :]=pd_raw.pd_rt[feeder,t+1:t+Pred_length].*pd_ratio;
    end

    da = pd_raw.pd_da[:,t+1+Pred_length:t+T];
    ratio = reshape(0.01:0.01/23:0.01/23*(Pred_length-1)+0.01, 1, Pred_length)
    pred_square = pred.^2;
    sigma_1 = hcat(ones(12,1)*ratio.*pred_square, 0.08*da.^2)
    traj_1 = hcat(pred, da);
    sigma_c=[];
    traj_c=[]
    for i=1:BN/3
        push!(sigma_c,sigma_1);
        push!(traj_c,traj_1);
    end
    traj=vcat(traj_c...)
    sigma=vcat(sigma_c...)
    # pd = pd_struct(traj, sigma);
    # pd = pd_struct(traj, sigma);
    pd = (traj=(traj), sigma = (sigma), ct=(ct))
    return pd
end

function pd_traj_pu_det(t, pd_raw, T, NoShunt, baseMVA)
    # this function generates the demand trajectory
    traj=zeros(NoShunt, 288)
    da = zeros(NoShunt, 288)
    qd_traj=zeros(NoShunt, 288)
    qd_da = zeros(NoShunt, 288)
    for shunt=1:NoShunt
        # println(multiplier.feeder_mult[bus])
        traj[shunt, :]=pd_raw.pd_rt[shunt, t+1:t+T]/baseMVA;
        da[shunt, :]=pd_raw.pd_da[shunt,t+1:t+T]/baseMVA;
        qd_traj[shunt, :]=pd_raw.qd_rt[shunt, t+1:t+T]/baseMVA;
        qd_da[shunt, :]=pd_raw.qd_da[shunt,t+1:t+T]/baseMVA;
    end
    ct=traj[:,t]
    qd_ct=qd_traj[:,t]
    sigma = zeros(NoShunt, T);
    # pd = pd_struct(traj, sigma);
    pd = (traj=(traj), sigma=(sigma), ct=(ct),
        qd_traj=(qd_traj), qd_ct=(qd_ct));
    return pd
end

function pg_traj(t, pg_raw, pg_noise, solar_error_max, p_rate, T, Pred_length);
    sg_max = maximum(p_rate*pg_raw.pg_rt, dims=2);
    # println(size(sg_max))
    rt_raw = p_rate*pg_raw.pg_rt;
    da_raw = p_rate*pg_raw.pg_da;
    mu_rt = rt_raw[:,t];
    mu_ct = rt_raw[:,t+1];
    pred = zeros(12, Pred_length)
    traj025 = 0.01:(0.025-0.01)/(Pred_length-1):0.025
    # traj = 0.01:(solar_error_max-0.01)/(Pred_length-1):solar_error_max;
    traj = traj025.+(solar_error_max-0.01)
    for feeder=1:12
        temp=sqrt.(traj./traj025);
        temp1=temp.*pg_noise[feeder][t,1:Pred_length];
        # pg_ratio = positive_array(pg_noise[feeder][t,1:Pred_length].+1);
        pg_ratio=positive_array(temp1.+1)
        pred[feeder, :]=rt_raw[feeder,t+1:t+Pred_length].*pg_ratio;
    end
    da = da_raw[:,t+1+Pred_length:t+T];
    mu_scenario = hcat(pred, da)
    mu = mu_scenario;
    # ratio = reshape(0.01:(solar_error_max-0.01)/23:(solar_error_max-0.01)/23*(Pred_length-1)+0.01, 1, Pred_length)
    ratio = reshape(traj025.+(solar_error_max-0.01),1,Pred_length)
    pred_square = pred.^2;
    sigma = hcat(ones(12,1)*ratio.*pred_square, p_rate^2*2*solar_error_max*da.^2);
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=(mu=(mu),sigma=(sigma), mu_ct=(mu_ct), sg_max=(sg_max));
    return pg
end

function pg_traj_pu_det(t, pg_raw, p_rate, T, NoShunt, baseMVA);
    #########
    #########
    sg_max = maximum(p_rate*pg_raw.pg_rt, dims=2);
    rt_raw = p_rate*pg_raw.pg_rt;
    da_raw = p_rate*pg_raw.pg_da;
    mu=zeros(NoShunt,T)
    da=zeros(NoShunt,T)
    for shunt = 1:NoShunt
        mu[shunt, :] = rt_raw[shunt,t+1:t+T]/baseMVA;
        da[shunt, :] = da_raw[shunt,t+1:t+T]/baseMVA;
    end

    mu_ct = mu[:,t+1];
    sigma = zeros(NoShunt,T)
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=(mu=(mu),mu_ct=(mu_ct),sigma=(sigma),sg_max=(sg_max));
    return pg
end




function read_price_data()
    filename = "../data/trace_feb28_2019.csv"
    data_trace = CSV.File(filename; dateformat="yyyy-mm-dd") |> DataFrame
    time_stamp = collect(data_trace[:,Symbol("RTD End Time Stamp")])
    ten_min_price = collect(data_trace[:,Symbol("RTD 10 Min Non Sync")])
    thirty_min_price = collect(data_trace[:,Symbol("RTD 30 Min Non Sync")])
    rt_lbmp = collect(data_trace[:, Symbol("Real Time LBMP")])
    da_lbmp = collect(data_trace[:, Symbol("Day Ahead LBMP")])
    price_raw = (timestamp = (time_stamp), RSRV_10=(ten_min_price),
        RSRV_30=(thirty_min_price), LMP_RT=(rt_lbmp), LMP_DA=(da_lbmp));
    return price_raw
end

###############################
function read_NY_demand_data(shunt)
    # P=shunt.P;
    NoShunt = length(shunt.P);
    demand_unit=matread("../data/NYISO-data/normalized_demand.mat")["normalized_demand"];
    #################################
    # 576*240
    # sum(demand_unit, dims=1) = ones
    #################################
    max_demand_unit = maximum(demand_unit[:,1:NoShunt], dims=1);

    frac_unit = reshape(ones(size(max_demand_unit))./max_demand_unit, NoShunt,1)
    frac = (shunt.P.*frac_unit*ones(1, 576))';
    demand_shunt = frac.*demand_unit[:, 1:NoShunt]
    # println(frac[1,:])

    frac_Q = (shunt.Q.*frac_unit*ones(1, 576))';
    reactive_demand_shunt = frac_Q.*demand_unit[:, 1:NoShunt]

    demand_da_unit=matread("../data/NYISO-data/normalized_demand_da.mat")["normalized_demand_da"]
    # max_demand_da_unit = maximum(demand_da_unit[:,1:NoShunt], dims=1);
    # frac_da_unit = reshape(ones(size(max_demand_da_unit))./max_demand_da_unit, NoShunt,1)
    # frac_da = reshape(shunt.P.*frac_da_unit*ones(1, 576), 576, NoShunt)
    demand_da_shunt = frac.*demand_da_unit[:, 1:NoShunt]
    # demand_da_shunt=10000*demand_da_unit[:, 1:NoShunt]
    # frac_da_Q = reshape(shunt.Q.*frac_da_unit*ones(1, 576), 576, NoShunt)
    reactive_demand_da_shunt = frac_Q.*demand_da_unit[:, 1:NoShunt]


    pd_raw = (demand_unit=(demand_unit), pd_rt = (demand_shunt'), pd_da = (demand_da_shunt'),
        qd_rt = (reactive_demand_shunt'), qd_da = (reactive_demand_da_shunt'));
end

###################################
# raw data
function read_NY_solar_data(pd_raw)
    # P=shunt.P;
    NoShunt = length(pd_raw.pd_rt[:,1]);
    solar_unit=matread("../data/NYISO-data/normalized_solar.mat")["normalized_solar"]
    #################################
    # 576*240
    # sum(demand_unit, dims=1) = ones
    #################################
    frac = ones(576,1)*
        (sum(pd_raw.pd_rt[:,1:288], dims=2)'
        ./sum(solar_unit[1:288, 1:NoShunt], dims=1));
    # println(size(
    # sum(pd_raw.pd_rt[:,1:288], dims=2)'./sum(solar_unit[1:288, 1:NoShunt], dims=1)
    # ))
    solar_shunt = frac.*solar_unit[:, 1:NoShunt];
    #
    solar_da_unit=
        matread("../data/NYISO-data/normalized_solar_da.mat")["normalized_solar_da"]
    solar_da_shunt = frac.*solar_da_unit[:, 1:NoShunt];

    pg_raw = (pg_rt = (solar_shunt'), pg_da = (solar_da_shunt'));
    return pg_raw
end
