# This runs all of the jobs necessary for our trade-off analysis

include("GML.jl")
using Plots
# # of timeslots
T=288;
# # of banks
BN=3;
# # of feeders
F=BN*4;
# # of secnarios
SN=6;

################################################################################

# Ancillar Markets Considered: "without", "10min", "30min"
ancillary_type = "10min"
# penetation levels
# p_rate = [0.25 0.5 0.75, 1]; #[0.25, 0.5, 0.75, 1];
p_rate = 0.75;

# [90.95,99]
# icdf = [-1.2816, -1.6449, -2.3263]
icdf = -1.2816;
# prediction length [2hr, 1hr, .5hr]
pred_length = 24;

# max solar error
solar_error_max_reservoir = [0.025 0.1 0.2];

# max solar error
B_cap = 3;
################################################################################

# load price data
price_raw = read_price_data()
delta_rt_raw=matread("../data/price_prediction.mat");

# read demand data
pd_raw = read_demand_data()
pd_noise = matread("../data/demand_noise.mat")["demand_noise"];

# read pg noise
pg_noise = matread("../data/solar_noise_0025.mat")["solar_noise"];

# read pg_raw
pg_raw = read_solar_data()

sigma_pd = zeros(288, length(solar_error_max_reservoir)*length(p_rate))
sigma_pg = zeros(288, length(solar_error_max_reservoir)*length(p_rate))
i=1;
for solar_error_max in solar_error_max_reservoir
    for current_time=1:T

        if current_time == 1
            global P_rsrv_feedback = [];
            global B_feedback=reshape([0.0910,0.0942,0.1434,0.1402,0.1382,0.2144,0.1192,0.1813,
                0.1355,0.0664,0.1182,0.0580],12,1);
        end

        feedback = (B_feedback=(B_feedback), P_rsrv_feedback=(P_rsrv_feedback));
        price = price_traj(current_time, ancillary_type, price_raw, delta_rt_raw, T, pred_length);
        # println(size(price.alpha_scenario))
        # println(size(price.lambda_scenario))
        pd = pd_traj(current_time, pd_raw, pd_noise, BN, T, pred_length);
        sigma_pd[current_time, i]=sum(pd.sigma[:,1]);
        # println(size(pd.traj))
        # println(size(pd.sigma))
        pg = pg_traj(current_time, pg_raw, pg_noise, solar_error_max, p_rate, T, pred_length);
        sigma_pg[current_time, i]=sum(pg.sigma[:,1]);
        # println(size(pg.mu))
        # println(size(pg.sigma))
        # obj = GML_Sys_Ava(T, F, BN, SN, pd, ancillary_type, icdf, B_cap);
        # val_opt = optimal_stoach_scenario(current_time, obj, feedback, pd, pg,
        #     price, ancillary_type);
        # for feeder = 1:12
        #     B_feedback[feeder, 1] = val_opt.B[feeder,1] - val_opt.R[feeder,1]/12;
        # end
        # if current_time == 1
        #     P_rsrv_feedback = [val_opt.P_rsrv]
        # else
        #     push!(P_rsrv_feedback,val_opt.P_rsrv)
        # end
        # mkpath(folder)
        # write_output_out(val_opt,
        #     string(folder, filename, "_time", current_time, ".csv"))
        # println("=================================================")
    end
    global i=i+1;
end
plot(1:288, sigma_pd[:,1], label="demand", linewidth=2)
plot!(1:288, sigma_pg[:,1], label="solar (0.025)", linewidth=2)

plot!(1:288, sigma_pg[:,2], label="solar (0.1)", linewidth=2)

plot!(1:288, sigma_pg[:,3], label="solar (0.2)", linewidth=2)
xlabel!("time")
ylabel!("variance")
# title!("New approach")
title!("Simulation 2")
