using CSV
using DataFrames
using Plots
using MAT


include("traj_gen_det.jl")

price_raw = read_price_data()
delta_rt_raw=matread("../data/price_prediction.mat");

## read read_demand_data
pd_raw = read_demand_data()
pd_noise = matread("../data/demand_noise.mat")["demand_noise"];
pg_raw = read_solar_data()

T = 288;

# Ancillar Markets Considered: "without", "10min", "30min"
a_types="10min"

# penetation levels
p_rate = 1

# chance constraint [50, 99] #[95, 90, 85]
cc = 90;

# prediction length
Pred_length = [20];

# max solar error
solar_error_max = [0.025, 0.05, 0.1, 0.5];

pd = pd_traj_det(1, pd_raw, T)
Pd_headnode = reshape(sum(pd.traj[1:4,:], dims=1),288,1);
pg = pg_traj_det(1, pg_raw, p_rate, T);
Pg_headnode_det1 = reshape(sum(pg.mu[1:4,:], dims=1),288,1);
Pg_headnode_det2 = reshape(sum(pg.mu[5:8,:], dims=1),288,1);

for p in p_rate
        folder = string("../../Energise-Tradeoff-git/GML_JHU_v2.2.0/results_0912/");
        folder = string(folder, "ten_min_anc/");
        filename = "ten_min_anc_";
        folder = string(folder, "p_rate/p_rate_", Integer(p*100), "/");
        filename = string(filename, "p_rate_", Integer(p*100), "_");
        folder = string(folder, "cc/cc_", cc, "/")
        filename = string(filename, "cc_", cc, "_")

        Cost = zeros(1,288)
        Pg = zeros(12,288)
        R = zeros(12,288)
        Qf = zeros(12,288)
        for time=1:T
            name=string(filename, "time", time, ".csv");
            data_trace = CSV.File(string(folder, name)) |> DataFrame
            Cost_temp= collect(data_trace[:,Symbol("Cost")])
            Cost[1,time] = Cost_temp[1]
            Pg_temp= collect(data_trace[:,Symbol("Pg")])
            Pg[:,time] = Pg_temp
            Qf_temp= collect(data_trace[:,Symbol("QF")])
            Qf[:,time] = Qf_temp
            R_temp= collect(data_trace[:,Symbol("R")])
            R[:,time] = R_temp
        end
        global Pg_headnode = reshape(sum(Pg[1:4,:], dims=1),288,1)
        global Qf_headnode = reshape(sum(Qf[1:4,:], dims=1),288,1)
        global R_headnode = reshape(sum(R[1:12,:], dims=1),288,1)
        global Pg_bank1 = reshape(sum(Pg[1:4,:], dims=1),288,1)
        global Pg_bank2 = reshape(sum(Pg[5:8,:], dims=1),288,1)
        global Pg_bank3 = reshape(sum(Pg[9:12,:], dims=1),288,1)
end
title_text = string("Chance Constraint Confidence ",cc, "%")
P=Pd_headnode-Pg_headnode-R_headnode;
Q=Qf_headnode;
plot(1:288, Pg_headnode, label="Pg bank1", lw=2)
plot!(1:288, Pg_bank2, label="Pg bank2", lw=2)
# plot!(1:288, sqrt(35*35/2)*ones(288,1), label="Pg Aval", lw=2)
# plot(1:288, Pd_headnode+Pg_headnode, label="Pg Aval", lw=2)
plot!(1:288, Pg_headnode_det1, label="Pg aval 1", lw=2)
plot!(1:288, Pg_headnode_det2, label="Pg aval 2", lw=2)
# plot!(1:288, Pd_headnode, label="Pd", lw=2)
# plot!(1:288, R_headnode, label="R", lw=2)
# plot!(1:288, Qf_headnode, label="Qf", lw=2)
# plot!(1:288, Pg_headnode_det-Pg_headnode, label="HN Cul", lw=2)
# plot!(1:288, Pg_bank1, label="B1", lw=2)
# plot!(1:288, Pg_bank2, label="B2", lw=2)
# plot!(1:288, Pg_bank3, label="B3", lw=2, title=title_text)
xlabel!("Time")
# ylabel!("Solar [MW]")

        # for pred_length in Pred_length
        #     folder = "results/"; filename = "";
        #     if a_type == "without"
        #         folder = string(folder, "no_anc/");
        #         filename = string(filename, "no_anc_");
        #     elseif a_type == "10min"
        #         folder = string(folder, "ten_min_anc/");
        #         filename = string(filename, "ten_min_anc_");
        #     end
        #     folder = string(folder, "p_rate/p_rate_", Integer(p*100), "/");
        #     filename = string(filename, "p_rate_", Integer(p*100), "_");
        #
        #     folder = string(folder, "pred/pred_", pred_length, "/")
        #     filename = string(filename, "pred_", pred_length, "_")
        #
        #     Cost = zeros(1,288)
        #     Pg = zeros(12,288)
        #     for time=1:T
        #         name=string(filename, "time", time, ".csv");
        #         data_trace = CSV.File(string(folder, name)) |> DataFrame
        #         Cost_temp= collect(data_trace[:,Symbol("Cost")])
        #         Cost[1,time] = Cost_temp[1]
        #         Pg_temp= collect(data_trace[:,Symbol("Pg")])
        #         Pg[:,time] = Pg_temp
        #     end
        #
        #     println(filename)
        #     println(sum(Pg))
        #     println(sum(Cost))
        #     println()
        # end

        # for sem in solar_error_max
        #
        #     folder = "results/"; filename = "";
        #     if a_type == "without"
        #         folder = string(folder, "no_anc/");
        #         filename = string(filename, "no_anc_");
        #     elseif a_type == "10min"
        #         folder = string(folder, "ten_min_anc/");
        #         filename = string(filename, "ten_min_anc_");
        #     end
        #     folder = string(folder, "p_rate/p_rate_", Integer(p*100), "/");
        #     filename = string(filename, "p_rate_", Integer(p*100), "_");
        #
        #     folder = string(folder, "solar/solar_", Integer(sem*1000), "/")
        #     filename = string(filename, "solar_", Integer(sem*1000), "_")
        #
        #     Cost = zeros(1,288)
        #     Pg = zeros(12,288)
        #     for time=1:T
        #         name=string(filename, "time", time, ".csv");
        #         data_trace = CSV.File(string(folder, name)) |> DataFrame
        #         Cost_temp= collect(data_trace[:,Symbol("Cost")])
        #         Cost[1,time] = Cost_temp[1]
        #         Pg_temp= collect(data_trace[:,Symbol("Pg")])
        #         Pg[:,time] = Pg_temp
        #         # Time_temp= collect(data_trace[:,Symbol("Time")])
        #         # Cost[1,time] = Time_temp[1]
        #     end
        #     println(filename)
        #     println(sum(Pg))
        #     println(sum(Cost))
        #     println()
        #
        # end
