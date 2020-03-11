using CSV
using DataFrames

T = 288;

# penetation levels
p_rate = [0.25, 0.5, 0.75, 1.0];

# chance constraint [50, 99] #[95, 90, 85]
cc = [90 95 99];

# prediction length
Pred_length = [30 60 120];

# max solar error
solar_error_max = [0.01, 0.02, 0.03, 0.035, 0.04, 0.05, 0.085];

B_cap = [3 15 30];
#
df_cost = DataFrame(CC=Integer[], pen_25=Float64[],
     pen_50=Float64[], pen_75=Float64[],pen_100=Float64[])
df_power = DataFrame(CC=Integer[], pen_25=Float64[],
     pen_50=Float64[], pen_75=Float64[],pen_100=Float64[])
for c in cc

    cost_sum_ra = [];
    pg_sum_ra = [];
    for p in p_rate

        folder = "results_1031/"; filename = "";
        folder = string(folder, "ten_min_anc/");
        filename = string(filename, "ten_min_anc_");
        folder = string(folder, "p_rate/p_rate_", Integer(p*100), "/");
        filename = string(filename, "p_rate_", Integer(p*100), "_");
        folder = string(folder, "cc/cc_", c, "/")
        filename = string(filename, "cc_", c, "_")

        Cost = zeros(1,288)
        Pg = zeros(12,288)
        for time=1:T
            name=string(filename, "time", time, ".csv");
            data_trace = CSV.File(string(folder, name)) |> DataFrame
            Cost_temp= collect(data_trace[:,Symbol("Cost")])
            Cost[1,time] = Cost_temp[1]
            Pg_temp= collect(data_trace[:,Symbol("Pg")])
            Pg[:,time] = Pg_temp
        end

        push!(cost_sum_ra, sum(Cost))
        push!(pg_sum_ra, sum(Pg))
    end
    push!(df_cost, [c; cost_sum_ra])
    push!(df_power, [c; pg_sum_ra])

end
mkpath("./results_1031/analysis/")
CSV.write("./results_1031/analysis/cc_table_cost.csv", df_cost)
CSV.write("./results_1031/analysis/cc_table_power.csv", df_power)

# df_cost = DataFrame(PredLength=Integer[], pen_25=Float64[],
#     pen_50=Float64[], pen_75=Float64[], pen_100=Float64[])
# df_power = DataFrame(PredLength=Integer[], pen_25=Float64[],
#     pen_50=Float64[], pen_75=Float64[], pen_100=Float64[])
# # for pred_length in Pred_length
# #     cost_sum_ra = [];
# #     pg_sum_ra = [];
# #     for p in p_rate
# #         folder = "./results_1011/"; filename = "";
# #         folder = string(folder, "ten_min_anc/");
# #         filename = string(filename, "ten_min_anc_");
# #         folder = string(folder, "p_rate/p_rate_", Integer(p*100), "/");
# #         filename = string(filename, "p_rate_", Integer(p*100), "_");
# #         folder = string(folder, "pred/pred_", pred_length, "/")
# #         filename = string(filename, "pred_", pred_length, "_")
# #
# #         Cost = zeros(1,288)
# #         Pg = zeros(12,288)
# #         for time=1:T
# #             name=string(filename, "time", time, ".csv");
# #             data_trace = CSV.File(string(folder, name)) |> DataFrame
# #             Cost_temp= collect(data_trace[:,Symbol("Cost")])
# #             Cost[1,time] = Cost_temp[1]
# #             Pg_temp= collect(data_trace[:,Symbol("Pg")])
# #             Pg[:,time] = Pg_temp
# #         end
# #         push!(cost_sum_ra, sum(Cost))
# #         push!(pg_sum_ra, sum(Pg))
# #     end
# #     push!(df_cost, [pred_length; cost_sum_ra])
# #     push!(df_power, [pred_length; pg_sum_ra])
# # end
# # CSV.write("./results_1011/analysis/pred_table_cost.csv", df_cost)
# # CSV.write("./results_1011/analysis/pred_table_power.csv", df_power)

df_cost = DataFrame(SEM=Float64[], pen_25=Float64[],
     pen_50=Float64[], pen_75=Float64[], pen_100=Float64[])
df_power = DataFrame(SEM=Float64[], pen_25=Float64[],
     pen_50=Float64[], pen_75=Float64[], pen_100=Float64[])
for sem in solar_error_max
    cost_sum_ra = [];
    pg_sum_ra = [];
    for p in p_rate

        folder = "./results_1031/"; filename = "";
        folder = string(folder, "ten_min_anc/");
        filename = string(filename, "ten_min_anc_");
        folder = string(folder, "p_rate/p_rate_", Integer(p*100), "/");
        filename = string(filename, "p_rate_", Integer(p*100), "_");
        folder = string(folder, "solar/solar_", Integer(sem*1000), "/")
        filename = string(filename, "solar_", Integer(sem*1000), "_")

        Cost = zeros(1,288)
        Pg = zeros(12,288)
        for time=1:T
            name=string(filename, "time", time, ".csv");
            data_trace = CSV.File(string(folder, name)) |> DataFrame
            Cost_temp= collect(data_trace[:,Symbol("Cost")])
            Cost[1,time] = Cost_temp[1]
            Pg_temp= collect(data_trace[:,Symbol("Pg")])
            Pg[:,time] = Pg_temp
            # Time_temp= collect(data_trace[:,Symbol("Time")])
            # Cost[1,time] = Time_temp[1]
        end
        push!(cost_sum_ra, sum(Cost))
        push!(pg_sum_ra, sum(Pg))
    end
    push!(df_cost, [sem; cost_sum_ra])
    push!(df_power, [sem; pg_sum_ra])
end
CSV.write("./results_1031/analysis/sem_table_cost.csv", df_cost)
CSV.write("./results_1031/analysis/sem_table_power.csv", df_power)


# df_cost = DataFrame(BCap=Integer[], pen_25=Float64[],
#     pen_50=Float64[], pen_75=Float64[], pen_100=Float64[])
# df_power = DataFrame(BCap=Integer[], pen_25=Float64[],
#     pen_50=Float64[], pen_75=Float64[], pen_100=Float64[])
# for b_cap in B_cap
#     cost_sum_ra = [];
#     pg_sum_ra = [];
#     for p in p_rate
#
#         folder = "./results_1011/"; filename = "";
#         folder = string(folder, "ten_min_anc/");
#         filename = string(filename, "ten_min_anc_");
#         folder = string(folder, "p_rate/p_rate_", Integer(p*100), "/");
#         filename = string(filename, "p_rate_", Integer(p*100), "_");
#         folder = string(folder, "b_cap/b_cap_", b_cap, "/")
#         filename = string(filename, "b_cap_", b_cap, "_")
#
#         Cost = zeros(1,288)
#         Pg = zeros(12,288)
#         for time=1:T
#             name=string(filename, "time", time, ".csv");
#             data_trace = CSV.File(string(folder, name)) |> DataFrame
#             Cost_temp= collect(data_trace[:,Symbol("Cost")])
#             Cost[1,time] = Cost_temp[1]
#             Pg_temp= collect(data_trace[:,Symbol("Pg")])
#             Pg[:,time] = Pg_temp
#             # Time_temp= collect(data_trace[:,Symbol("Time")])
#             # Cost[1,time] = Time_temp[1]
#         end
#         push!(cost_sum_ra, sum(Cost))
#         push!(pg_sum_ra, sum(Pg))
#     end
#     push!(df_cost, [b_cap; cost_sum_ra])
#     push!(df_power, [b_cap; pg_sum_ra])
# end
# CSV.write("./results_1011/analysis/bat_table_cost.csv", df_cost)
# CSV.write("./results_1011/analysis/bat_table_power.csv", df_power)
return "finish"
