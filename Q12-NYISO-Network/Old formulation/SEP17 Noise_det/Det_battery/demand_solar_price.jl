function price_traj(t, ancillary_type, price_raw, delta_rt_raw, T, Pred_length)
# this function generates the price trajectory
    PTID=3;
    zero_flag = 1;
    Pred_scen=3;

    if ancillary_type == "without"
        alpha_ct=0;
        alpha_scenario=zeros(Pred_scen,T);
    elseif ancillary_type == "10min"
        alpha_rt=price_raw.RSRV_10[t];
        alpha_ct=price_raw.RSRV_10[t+1];
        if alpha_rt == 0
            zero_flag = 1
        else
            zero_flag = 2
        end
        alpha_pd=
            delta_rt_raw["RSRV10centroid"][PTID,t,zero_flag][:,1:Pred_length];
        alpha_scenario=hcat(alpha_pd,zeros(Pred_scen,T-Pred_length));
    elseif ancillary_type == "30min"
        alpha_rt=price_raw.RSRV_30[t];
        alpha_ct=price_raw.RSRV_30[t+1];
        if alpha_rt == 0
            zero_flag = 1
        else
            zero_flag = 2
        end
        alpha_pd=
            delta_rt_raw["RSRV30centroid"][PTID,t,zero_flag][:,1:Pred_length];
        alpha_scenario=hcat(alpha_pd,zeros(Pred_scen,T-Pred_length));
    end

    # println(zero_flag)

    lambda_rt=price_raw.LMP_RT[t]
    lambda_ct=price_raw.LMP_RT[t+1]
    # lambda_rt=price_raw["lambda_rt"][t]
    lambda_pd=lambda_rt.+delta_rt_raw["RTcentroid"][PTID,t,zero_flag][:,1:Pred_length];
    # Pred_scen=length(lambda_pd[:,1]);

    # println(size(lambda_pd))
    lambda_da=price_raw.LMP_DA[t+Pred_length+1:t+T];
    # lambda_da=price_raw["lambda_da"][t+Pred_length+1:t+T-1];
    # println(size(lambda_da))
    lambda_offline=ones(Pred_scen,1)*reshape(lambda_da, 1, length(lambda_da));
    lambda_scenario = hcat(lambda_pd,lambda_offline);



    probability=delta_rt_raw["probability"][PTID,t,zero_flag];
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
    lambda_scenario = ones(3,1)*reshape(price_raw.LMP_RT[t+1:t+T], 1, T);
    if ancillary_type == "without"
        alpha_ct=0;
        alpha_scenario=zeros(3,T);
    elseif ancillary_type == "10min"
        alpha_ct=price_raw.RSRV_10[t+1];
        # println(alpha_ct)
        alpha_scenario=ones(3,1)*reshape(price_raw.RSRV_10[t+1:t+T], 1, T);
    elseif ancillary_type == "30min"
        alpha_ct=price_raw.RSRV_30[t+1];
        alpha_scenario=ones(3,1)*reshape(price_raw.RSRV_30[t+1:t+T], 1, T);
    end
    probability=[1 0 0];
    price = (lambda_ct=(lambda_ct), lambda_scenario = (lambda_scenario),
        alpha_ct=(alpha_ct), alpha_scenario=(alpha_scenario),
        probability=(probability))
    return price;
end

function pd_traj(date, t, pd_raw, pd_noise, T, NoShunt, Pred_length, baseMVA)

    feeder_real = ORU_feeder_Pd(date);
    feeder_reactive = ORU_feeder_Qd(date);

    traj=zeros(NoShunt, 288)
    da = zeros(NoShunt, 288-Pred_length)
    qd_traj=zeros(NoShunt, 288)
    qd_da = zeros(NoShunt, 288)

    raw_pu_traj=zeros(NoShunt, 288);
    raw_pu_da = zeros(NoShunt, 288);
    raw_pu_qd_traj=zeros(NoShunt, 288);
    raw_pu_qd_da = zeros(NoShunt, 288);

    for shunt=1:NoShunt
        # println(multiplier.feeder_mult[bus])
        if shunt == 2
            raw_pu_traj[shunt, :]=feeder_real.feeder_pd[2, t+1:t+T]/baseMVA;
            raw_pu_qd_traj[shunt, :]=feeder_reactive.feeder_qd[2, t+1:t+T]/baseMVA;

            raw_pu_da[shunt, :]=feeder_real.feeder_pd_da[2, t+1:t+T]/baseMVA;
            raw_pu_qd_da[shunt, :]=feeder_reactive.feeder_qd_da[2,t+1:t+T]/baseMVA;
        elseif shunt == 4
            raw_pu_traj[shunt, :]=feeder_real.feeder_pd[1, t+1:t+T]/baseMVA;
            raw_pu_qd_traj[shunt, :]=feeder_reactive.feeder_qd[1, t+1:t+T]/baseMVA;

            raw_pu_da[shunt, :]=feeder_real.feeder_pd_da[1, t+1:t+T]/baseMVA;
            raw_pu_qd_da[shunt, :]=feeder_reactive.feeder_qd_da[1,t+1:t+T]/baseMVA;
        elseif shunt == 5
            raw_pu_traj[shunt, :]=feeder_real.feeder_pd[3, t+1:t+T]/baseMVA;
            raw_pu_qd_traj[shunt, :]=feeder_reactive.feeder_qd[3, t+1:t+T]/baseMVA;

            raw_pu_da[shunt, :]=feeder_real.feeder_pd_da[3, t+1:t+T]/baseMVA;
            raw_pu_qd_da[shunt, :]=feeder_reactive.feeder_qd_da[3,t+1:t+T]/baseMVA;
        else
            raw_pu_traj[shunt, :]=pd_raw.pd_rt[shunt, t+1:t+T]/baseMVA;
            raw_pu_qd_traj[shunt, :]=pd_raw.qd_rt[shunt, t+1:t+T]/baseMVA;
            raw_pu_da[shunt, :]=pd_raw.pd_da[shunt,t+1:t+T]/baseMVA;
            raw_pu_qd_da[shunt, :]=pd_raw.qd_da[shunt,t+1:t+T]/baseMVA;
        end

    end

    ct=raw_pu_traj[:, 1]
    qd_ct=raw_pu_qd_traj[:, 1]

    # pd = pd_struct(traj, sigma);
    # pd = (traj=(traj), sigma=(sigma), ct=(ct),
    #     qd_traj=(qd_traj), qd_ct=(qd_ct));
    # return pd

    pred = zeros(NoShunt, Pred_length);
    qd_pred = zeros(NoShunt, Pred_length);
    for shunt=1:NoShunt
        noise_flag = mod(shunt, 12)+1;
        pd_ratio = positive_array(pd_noise[noise_flag][t,1:Pred_length].+1);
        pred[shunt, :]=raw_pu_traj[shunt,1:Pred_length].*pd_ratio;
        qd_pred[shunt, :]=raw_pu_qd_traj[shunt,1:Pred_length].*pd_ratio;
    end

    ratio = reshape(0.01:0.01/23:0.01/23*(Pred_length-1)+0.01, 1, Pred_length)
    pred_square = pred.^2;
    qd_pred_square = qd_pred.^2;

    sigma = hcat(ones(NoShunt,1)*ratio.*pred_square, 0.08*da.^2)
    qd_sigma = hcat(ones(NoShunt,1)*ratio.*qd_pred_square, 0.08*da.^2)

    da=raw_pu_da[:,1+Pred_length:T]
    qd_da = raw_pu_qd_da[:,1+Pred_length:T]

    traj = hcat(pred, da);
    qd_traj = hcat(qd_pred, qd_da);

    pd = (traj=(traj), sigma=(sigma), ct=(ct),
        qd_traj=(qd_traj), qd_sigma=(qd_sigma), qd_ct=(qd_ct));
    return pd
end

function pd_traj_pu_det(date, t, pd_raw, T, NoShunt, baseMVA)
    # this function generates the demand trajectory
    traj=zeros(NoShunt, 288)
    da = zeros(NoShunt, 288)
    qd_traj=zeros(NoShunt, 288)
    qd_da = zeros(NoShunt, 288)

    feeder_real = ORU_feeder_Pd(date);
    feeder_reactive = ORU_feeder_Qd(date);

    raw_pu_traj=zeros(NoShunt, 288);
    raw_pu_da = zeros(NoShunt, 288);
    raw_pu_qd_traj=zeros(NoShunt, 288);
    raw_pu_qd_da = zeros(NoShunt, 288);

    for shunt=1:NoShunt
        # println(multiplier.feeder_mult[bus])
        if shunt == 2
            raw_pu_traj[shunt, :]=feeder_real.feeder_pd[2, t+1:t+T]/baseMVA;
            raw_pu_qd_traj[shunt, :]=feeder_reactive.feeder_qd[2, t+1:t+T]/baseMVA;

            raw_pu_da[shunt, :]=feeder_real.feeder_pd_da[2, t+1:t+T]/baseMVA;
            raw_pu_qd_da[shunt, :]=feeder_reactive.feeder_qd_da[2,t+1:t+T]/baseMVA;
        elseif shunt == 4
            raw_pu_traj[shunt, :]=feeder_real.feeder_pd[1, t+1:t+T]/baseMVA;
            raw_pu_qd_traj[shunt, :]=feeder_reactive.feeder_qd[1, t+1:t+T]/baseMVA;

            raw_pu_da[shunt, :]=feeder_real.feeder_pd_da[1, t+1:t+T]/baseMVA;
            raw_pu_qd_da[shunt, :]=feeder_reactive.feeder_qd_da[1,t+1:t+T]/baseMVA;
        elseif shunt == 5
            raw_pu_traj[shunt, :]=feeder_real.feeder_pd[3, t+1:t+T]/baseMVA;
            raw_pu_qd_traj[shunt, :]=feeder_reactive.feeder_qd[3, t+1:t+T]/baseMVA;

            raw_pu_da[shunt, :]=feeder_real.feeder_pd_da[3, t+1:t+T]/baseMVA;
            raw_pu_qd_da[shunt, :]=feeder_reactive.feeder_qd_da[3,t+1:t+T]/baseMVA;
        else
            raw_pu_traj[shunt, :]=pd_raw.pd_rt[shunt, t+1:t+T]/baseMVA;
            raw_pu_qd_traj[shunt, :]=pd_raw.qd_rt[shunt, t+1:t+T]/baseMVA;
            raw_pu_da[shunt, :]=pd_raw.pd_da[shunt,t+1:t+T]/baseMVA;
            raw_pu_qd_da[shunt, :]=pd_raw.qd_da[shunt,t+1:t+T]/baseMVA;
        end

    end

    for shunt=1:NoShunt
        # println(multiplier.feeder_mult[bus])
        traj[shunt, :]=raw_pu_traj[shunt, :];
        da[shunt, :]=raw_pu_da[shunt,:];
        qd_traj[shunt, :]=raw_pu_qd_traj[shunt, :];
        qd_da[shunt, :]=raw_pu_qd_da[shunt,:];
    end
    ct=traj[:,1]
    qd_ct=qd_traj[:,1]
    sigma = zeros(NoShunt, T);
    qd_sigma = zeros(NoShunt, T);
    # pd = pd_struct(traj, sigma);
    pd = (traj=(traj), sigma=(sigma), ct=(ct),
        qd_traj=(qd_traj), qd_sigma = (qd_sigma), qd_ct=(qd_ct));
    return pd
end

function pg_traj(date, t, pg_raw, pg_noise, solar_error_max, p_rate, T, Pred_length,
    NoShunt, baseMVA);
    feeder_solar = ORU_feeder_Pg(date);
    sg_max = maximum(p_rate*pg_raw.pg_rt, dims=2)/baseMVA;

    sg_max[2] = 5.4/baseMVA;
    sg_max[4] = 34.655/baseMVA;
    sg_max[5] = 13.225/baseMVA;

    rt_raw = p_rate*pg_raw.pg_rt/baseMVA;
    rt_raw[2,:] = feeder_solar.feeder_pg[2,:]/baseMVA;
    rt_raw[4,:] = feeder_solar.feeder_pg[1,:]/baseMVA;
    rt_raw[5,:] = feeder_solar.feeder_pg[3,:]/baseMVA;

    da_raw = p_rate*pg_raw.pg_da/baseMVA;
    da_raw[2,:] = feeder_solar.feeder_pg_da[2,:]/baseMVA;
    da_raw[4,:] = feeder_solar.feeder_pg_da[1,:]/baseMVA;
    da_raw[5,:] = feeder_solar.feeder_pg_da[3,:]/baseMVA;

    # sg_max[2] = 5.4;
    # sg_max[4] = 34.655;
    # sg_max[5] = 13.225;

    # rt_raw = p_rate*pg_raw.pg_rt/baseMVA;
    # rt_raw[2,:] = feeder_solar.feeder_pg[2,:];
    # rt_raw[4,:] = feeder_solar.feeder_pg[1,:];
    # rt_raw[5,:] = feeder_solar.feeder_pg[3,:];

    # da_raw = p_rate*pg_raw.pg_da/baseMVA;
    # da_raw[2,:] = feeder_solar.feeder_pg_da[2,:];
    # da_raw[4,:] = feeder_solar.feeder_pg_da[1,:];
    # da_raw[5,:] = feeder_solar.feeder_pg_da[3,:];

    mu_rt = rt_raw[:,t];
    mu_ct = rt_raw[:,t+1];

    if solar_error_max >=0

        pred = zeros(NoShunt, Pred_length);
        qd_pred = zeros(NoShunt, Pred_length);

        traj025 = 0.01:(0.025-0.01)/(Pred_length-1):0.025
        # traj = 0.01:(solar_error_max-0.01)/(Pred_length-1):solar_error_max;
        traj = traj025.+(solar_error_max-0.01)

        for shunt=1:NoShunt
            noise_flag = mod(shunt, 12)+1;

            solar_error_mult=sqrt.(traj./traj025);
            solar_error_std=solar_error_mult.*pg_noise[noise_flag][t,1:Pred_length];
            # pg_ratio = positive_array(pg_noise[feeder][t,1:Pred_length].+1);
            pg_ratio=positive_array(solar_error_std.+1)
            pred[shunt, :]=rt_raw[shunt,t+1:t+Pred_length].*pg_ratio;
        end

        da = da_raw[:,t+1+Pred_length:t+T];
        # ratio = reshape(0.01:(solar_error_max-0.01)/23:(solar_error_max-0.01)/23*(Pred_length-1)+0.01, 1, Pred_length)
        ratio = reshape(traj025.+(solar_error_max-0.01),1,Pred_length)
        pred_square = pred.^2;
        sigma = hcat(ones(NoShunt,1)*ratio.*pred_square, p_rate^2*2*solar_error_max*da.^2);
    else
        pred = rt_raw[:, t:t+Pred_length];
        da = da_raw[:,t+1+Pred_length:t+T];
        sigma = zeros(NoShunt, T);
    end

    mu = hcat(pred, da)
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=(mu=(mu),sigma=(sigma), mu_ct=(mu_ct), sg_max=(sg_max));
    return pg
end

function pg_traj_pu_det(date, t, pg_raw, p_rate, T, NoShunt, baseMVA);
    #########
    #########
    feeder_solar = ORU_feeder_Pg(date);
    sg_max = maximum(p_rate*pg_raw.pg_rt, dims=2)/baseMVA;

    sg_max[2] = 5.4/baseMVA;
    sg_max[4] = 34.655/baseMVA;
    sg_max[5] = 13.225/baseMVA;

    rt_raw = p_rate*pg_raw.pg_rt/baseMVA;
    rt_raw[2,:] = feeder_solar.feeder_pg[2,:]/baseMVA;
    rt_raw[4,:] = feeder_solar.feeder_pg[1,:]/baseMVA;
    rt_raw[5,:] = feeder_solar.feeder_pg[3,:]/baseMVA;

    da_raw = p_rate*pg_raw.pg_da/baseMVA;
    da_raw[2,:] = feeder_solar.feeder_pg_da[2,:]/baseMVA;
    da_raw[4,:] = feeder_solar.feeder_pg_da[1,:]/baseMVA;
    da_raw[5,:] = feeder_solar.feeder_pg_da[3,:]/baseMVA;

    # sg_max[2] = 5.4;
    # sg_max[4] = 34.655;
    # sg_max[5] = 13.225;

    # rt_raw = p_rate*pg_raw.pg_rt/baseMVA;
    # rt_raw[2,:] = feeder_solar.feeder_pg[2,:];
    # rt_raw[4,:] = feeder_solar.feeder_pg[1,:];
    # rt_raw[5,:] = feeder_solar.feeder_pg[3,:];

    # da_raw = p_rate*pg_raw.pg_da/baseMVA;
    # da_raw[2,:] = feeder_solar.feeder_pg_da[2,:];
    # da_raw[4,:] = feeder_solar.feeder_pg_da[1,:];
    # da_raw[5,:] = feeder_solar.feeder_pg_da[3,:];

    mu=zeros(NoShunt,T)
    da=zeros(NoShunt,T)
    for shunt = 1:NoShunt
        mu[shunt, :] = rt_raw[shunt,t+1:t+T];
        da[shunt, :] = da_raw[shunt,t+1:t+T];
    end

    mu_ct = mu[:,1];
    sigma = zeros(NoShunt,T)
    # println(sum(sg_max))
    # println(sum(mu))
    # pg=pg_struct(mu,mu_rt,mu_scenario,sigma);
    pg=(mu=(mu),mu_ct=(mu_ct),sigma=(sigma),sg_max=(sg_max));
    return pg
end




function read_price_data(date)
    filename = string("GML_data/price-solar-demand/aug_traj/trace_",date,"_2019.csv")
    # filename = "GML_data/price-solar-demand/trace_feb28_2019.csv"
    data_trace = CSV.File(filename; dateformat="yyyy-mm-dd") |> DataFrame
    time_stamp = collect(data_trace[:,Symbol("RTD_End_Time_Stamp")])
    ten_min_price = collect(data_trace[:,Symbol("RTD_10_Min_Non_Sync")])
    thirty_min_price = collect(data_trace[:,Symbol("RTD_30_Min_Non_Sync")])
    rt_lbmp = collect(data_trace[:, Symbol("Real_Time_LBMP")])
    da_lbmp = collect(data_trace[:, Symbol("Day_Ahead_LBMP")])
    price_raw = (timestamp = (time_stamp), RSRV_10=(ten_min_price),
        RSRV_30=(thirty_min_price), LMP_RT=(rt_lbmp), LMP_DA=(da_lbmp));
    return price_raw
end

###############################
function read_NY_demand_data(shunt)
    # P=shunt.P;
    NoShunt = length(shunt.P);
    demand_unit=matread("GML_data/price-solar-demand/normalized_demand.mat")["normalized_demand"];
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

    demand_da_unit=matread("GML_data/price-solar-demand/normalized_demand_da.mat")["normalized_demand_da"]
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
    solar_unit=matread("GML_data/price-solar-demand/normalized_solar.mat")["normalized_solar"]
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
        matread("GML_data/price-solar-demand/normalized_solar_da.mat")["normalized_solar_da"]
    solar_da_shunt = frac.*solar_da_unit[:, 1:NoShunt];

    pg_raw = (pg_rt = (solar_shunt'), pg_da = (solar_da_shunt'));
    return pg_raw
end

function pre_process(date)
    data_trace_bus=CSV.File("GML_data/NYISO-data/bus.csv") |> DataFrame
    data_trace_shunt=CSV.File("GML_data/NYISO-data/shunt_ORU.csv") |> DataFrame
    data_trace_branch=CSV.File("GML_data/NYISO-data/branch_edit.csv") |> DataFrame
    data_trace_gen=CSV.File("GML_data/NYISO-data/gen.csv") |> DataFrame



    ################################################################################
    shunt_struct = NYISO_shunt_data(data_trace_bus, data_trace_shunt);
    bus_struct = NYISO_bus_data(data_trace_bus);
    branch_struct =NYISO_branch_data(data_trace_branch);
    gen_struct = NYISO_gen_data(data_trace_gen);

    ################################################################################
    ## load raw demand and price data
    # # load price data
    # In file demand_solar_price
    price_raw = read_price_data(date)
    # delta_rt_raw=matread("../data/price_prediction.mat");
    delta_rt_raw=matread("GML_data/price-solar-demand/Scenarios.mat");
    # println(size(delta_rt_raw["probability"]))

    # ## read read_demand_data
    pd_raw = read_NY_demand_data(shunt_struct)
    # # println(size(pd_raw.pd_rt))
    # # plot(reshape(pd_raw.pd_rt[1,:],576,1),label="pd rt", linewidth=2)
    # # plot!(reshape(pd_raw.qd_rt[1,:],576,1),label="Qd rt", linewidth=2)
    # # println(sum(pd_raw.pd_rt, dims=1))
    # # println(sum(pd_raw.pd_rt, dims=2))
    # # println(sum(pd_raw.pd_rt))
    #
    #
    # ###########################################################################
    pd_noise = matread("GML_data/price-solar-demand/demand_noise.mat")["demand_noise"]/10;
    # println(size(pd_noise))
    # println(size(pd_noise[1]))

    #
    #
    #
    # # ###########################################################################
    pg_raw = read_NY_solar_data(pd_raw)
    # # #
    # # # println(size(pg_raw.pg_rt))
    # # # println(sum(pg_raw.pg_rt))
    # # # println(size(pg_raw.pg_da))
    # # # println(sum(pg_raw.pg_da))
    # # # plot(1:576, reshape(sum(pd_raw.pd_rt, dims=1),576,1), label="pd rt", linewidth=2)
    # # # plot!(1:576, reshape(sum(pd_raw.pd_da, dims=1),576,1), label="pd da", linewidth=2)
    # # # plot!(1:576, reshape(sum(pg_raw.pg_rt, dims=1),576,1), label="pg rt", linewidth=2)
    # # # plot!(1:576, reshape(sum(pg_raw.pg_da, dims=1),576,1), label="pg da", linewidth=2)
    # #
    # # ###########################################################################
    return raw_data = (
        shunt_struct = (shunt_struct),
        bus_struct = (bus_struct),
        branch_struct = (branch_struct),
        gen_struct = (gen_struct),
        price_raw = (price_raw),
        delta_rt_raw = (delta_rt_raw),
        pd_raw = (pd_raw),
        pd_noise  = (pd_noise),
        pg_raw = (pg_raw)
        )
end

function ORU_feeder_Pd(date)
    feeder_list = ["Feeder1","Feeder2","Feeder3"]
    feeder_pd = zeros(3, 576)
    for feeder = 1:3
        filename = string("GML_data/ORU_feeder/",feeder_list[feeder],"_P/",date,".csv");
        data_trace_feeder=CSV.File(filename) |> DataFrame
        feeder_pd[feeder, :] = reshape(data_trace_feeder.Pd_MW, 1, 576)
    end
    feeder_pd_da = zeros(3, 576)
    for time =1:576
        # hour_step = ceil(Int, time/12)
        hour_pos = (ceil(Int, time/12)-1)*12+1;
        feeder_pd_da[:, time] =feeder_pd[:, hour_pos];
    end
    return feeder_real_demand = (feeder_pd=(feeder_pd), feeder_pd_da = (feeder_pd_da));
end

function ORU_feeder_Pg(date)
    feeder_list = ["Feeder1","Feeder2","Feeder3"]
    feeder_pg = zeros(3, 576)
    for feeder = 1:3
        filename = string("GML_data/ORU_feeder/",feeder_list[feeder],"_Pg/",date,".csv");
        data_trace_feeder=CSV.File(filename) |> DataFrame
        feeder_pg[feeder, :] = reshape(data_trace_feeder.solar_MW, 1, 576)
    end
    feeder_pg_da = zeros(3, 576)
    for time =1:576
        # hour_step = ceil(Int, time/12)
        hour_pos = (ceil(Int, time/12)-1)*12+1;
        feeder_pg_da[:, time] =feeder_pg[:, hour_pos];
    end
    return feeder_real_solar = (feeder_pg=(feeder_pg), feeder_pg_da = (feeder_pg_da));
end

function ORU_feeder_Qd(date)
    feeder_list = ["Feeder1","Feeder2","Feeder3"]
    feeder_qd = zeros(3, 576)
    for feeder = 1:3
        filename = string("GML_data/ORU_feeder/",feeder_list[feeder],"_Q/",date,".csv");
        data_trace_feeder=CSV.File(filename) |> DataFrame
        feeder_qd[feeder, :] = reshape(data_trace_feeder.Qd_MVAR, 1, 576)
    end

    feeder_qd_da = zeros(3, 576)
    for time =1:576
        # hour_step = ceil(Int, time/12)
        hour_pos = (ceil(Int, time/12)-1)*12+1;
        feeder_qd_da[:, time] =feeder_qd[:, hour_pos];
    end
    return feeder_reactive_demand = (feeder_qd=(feeder_qd), feeder_qd_da = (feeder_qd_da));
end
