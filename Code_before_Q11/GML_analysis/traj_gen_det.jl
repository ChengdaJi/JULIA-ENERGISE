function price_traj_det(t, ancillary_type, price_raw, T)
# this function generates the price trajectory
    PTID=3;
    lambda_rt=price_raw.LMP_RT[t]
    lambda_scenario = ones(6,1)*reshape(price_raw.LMP_RT[t+1:t+T-1], 1, T-1);
    if ancillary_type == "without"
        alpha_rt=0;
        alpha_scenario=zeros(6,T-1);
    elseif ancillary_type == "10min"
        alpha_rt=price_raw.RSRV_10[t];
        alpha_scenario=ones(6,1)*reshape(price_raw.RSRV_10[t+1:t+T-1], 1, T-1);
    elseif ancillary_type == "30min"
        alpha_rt=price_raw.RSRV_30[t];
        alpha_scenario=ones(6,1)*reshape(price_raw.RSRV_30[t+1:t+T-1], 1, T-1);
    end
    probability=[1 0 0 0 0 0];
    price = price_struct(lambda_rt, lambda_scenario, alpha_rt, alpha_scenario, probability)
    return price;
end

function pd_traj_det(t, pd_raw, T)
    # this function generates the demand trajectory
    traj=pd_raw.pd_rt[:, t:t+T-1];
    da=pd_raw.pd_da[:,t:t+T-1];
    sigma = zeros(12, T-1);
    # pd = pd_struct(traj, sigma);
    pd = (da=(da), traj=(traj), sigma=(sigma));
    return pd
end

function pg_traj_det(t, pg_raw, p_rate, T);
    rt_raw = p_rate*pg_raw.pg_rt;
    da_raw = p_rate*pg_raw.pg_da;
    mu = rt_raw[:,t:t+T-1];
    da = da_raw[:,t:t+T-1];
    mu_rt = rt_raw[:,t];
    mu_scenario=rt_raw[:,t+1:t+T-1];
    sigma = zeros(12,T-1)
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=(da=(da),mu=(mu),mu_rt=(mu_rt),mu_scenario=(mu_scenario),sigma=(sigma));
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
        RSRV_30=(ten_min_price), LMP_RT=(rt_lbmp), LMP_DA=(da_lbmp));
    return price_raw
end

function read_demand_data()
    filename = "../data/demand_feb_6_7.csv"
    data_trace = CSV.File(filename; dateformat="yyyy-mm-dd") |> DataFrame
    pd_da_sum= collect(data_trace[:,Symbol("da_sum")])
    # pd_feeder_2 = collect(data_trace[:,Symbol("ckt2_MW_bg276")])
    # pd_rt = hcat(pd_feeder_1, pd_feeder_2)
    pd_rt=zeros(12, 576)
    pd_da=zeros(12, 576)
    for time=1:576
        for feeder=1:12
            pd_rt[feeder, time] = data_trace[time,feeder+2];
            pd_da[feeder, time] = data_trace[time,feeder+15];
        end
    end
    pd_raw = (pd_rt = (pd_rt), pd_da_sum =(pd_da_sum), pd_da = (pd_da));
    return pd_raw
end


function read_solar_data()
    filename = "../data/Solar.csv"
    data_trace = CSV.File(filename; dateformat="yyyy-mm-dd") |> DataFrame

    # pd_feeder_2 = collect(data_trace[:,Symbol("ckt2_MW_bg276")])
    # pd_rt = hcat(pd_feeder_1, pd_feeder_2)
    pg_rt=zeros(12, 576)
    pg_da=zeros(12, 576)
    for time=1:576
        for feeder=1:12
            pg_rt[feeder, time] = data_trace[time,feeder+1];
            pg_da[feeder, time] = data_trace[(ceil(Int8, time/12))*12,feeder+1];
        end
    end
    pg_raw = (pg_rt = (pg_rt), pg_da = (pg_da));
    return pg_raw
end
