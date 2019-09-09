using CSV
using DataFrames

Cost = zeros(1,288)
Pg = zeros(12,288)
for time=1:T
    name=string("results/Yue/pred_length/2hr/Time", time, ".csv");
    data_trace = CSV.File(name) |> DataFrame
    Cost_temp= collect(data_trace[:,Symbol("Cost")])
    Cost[1,time] = Cost_temp[1]
    Pg_temp= collect(data_trace[:,Symbol("Pg")])
    Pg[:,time] = Pg_temp
end
println(sum(Pg))
println(sum(Cost))
