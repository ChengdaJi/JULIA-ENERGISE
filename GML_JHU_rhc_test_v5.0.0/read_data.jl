using MAT
using CSV
using DataFrames


# data_trace=CSV.File("../data/NYISO-data/shunt.csv") |> DataFrame
# P = collect(data_trace[:,Symbol("P")]).+0.01;


# demand_unit=matread("../data/NYISO-data/normalized_demand.mat")["normalized_demand"]
# demand=zeros(576, length(P))
# for n_shunt=1:length(P)
#     demand[:, n_shunt] = demand_unit[:, n_shunt].*(576).*P[n_shunt].*0.7;
# end
# demand_time = sum(demand, dims=2);
# println(maximum(demand_time))
# println(sum(demand_time))
#
# demand_da_unit=matread("../data/NYISO-data/normalized_demand_da.mat")["normalized_demand_da"]
# demand_da=zeros(576, length(P))
# for n_shunt=1:length(P)
#     demand_da[:, n_shunt] = demand_da_unit[:, n_shunt].*(576).*P[n_shunt].*0.7;
# end
# demand_da_time = sum(demand_da, dims=2);
# println(maximum(demand_da_time))
# println(sum(demand_da_time))
#
# solar_unit=matread("../data/NYISO-data/normalized_solar.mat")["normalized_solar"]
# solar=zeros(576, length(P))
# for n_shunt=1:length(P)
#     solar[:, n_shunt] = solar_unit[:, n_shunt].*(576).*P[n_shunt].*0.7;
# end
# solar_time = sum(solar, dims=2);
# println(maximum(solar_time))
#
# solar_da_unit=matread("../data/NYISO-data/normalized_solar_da.mat")["normalized_solar_da"]
# solar_da=zeros(576, length(P))
# for n_shunt=1:length(P)
#     solar_da[:, n_shunt] = solar_da_unit[:, n_shunt].*(576).*P[n_shunt].*0.7;
# end
# solar_da_time = sum(solar_da, dims=2);
# println(maximum(solar_da_time))
function shunt_data()
    data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
    baseKV_bus = collect(data_trace_bus[:,Symbol("baseKV")]);
    println(baseKV_bus[64])
    data_trace_shunt=CSV.File("../data/NYISO-data/shunt.csv") |> DataFrame
    find_bus = collect(data_trace_shunt[:,Symbol("bus")]);
    Q_max = collect(data_trace_shunt[:,Symbol("Q_max")]);
    baseKV = zeros(209,1)
    for nos = 1:lenght(of_bus)
        baseKV[nos]=baseKV_bus[find_bus[nos]]
    end
    shunt = (find_bus = (find_bus), baseKV=(baseKV), Q_max=(Q_max))
end
