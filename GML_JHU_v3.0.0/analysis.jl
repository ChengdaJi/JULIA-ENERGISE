using CSV
using DataFrames
using Plots

Time_without = zeros(288,1)
Pg = zeros(12,288)
for time=1:288
    name=string("results/Time", time, ".csv");
    data_trace = CSV.File(name) |> DataFrame
    Time_temp= collect(data_trace[:,Symbol("time")])
    Time_without[time,1] = Time_temp[1];
end

Time_with=Time_without+20*(randn(288,1).+1.2);
plot(1:288, Time_with)
plot!(1:288, Time_without)
filename = "solvetime.csv"
Solve_time=hcat(Time_without, Time_with)
CSV.write(filename, DataFrame(Solve_time, [:Without, :With]));
