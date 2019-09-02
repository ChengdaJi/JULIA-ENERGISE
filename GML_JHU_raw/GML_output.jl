function GML_output(opt_value)
###############################################################################
# write the solar file
println("GML - Start writing files!");
total_n=size(opt_value.Pg)[2];
time_stamp=[1:total_n]
# Vector(Int64, size(opt_value.Pg)[2])
time_stamp=zeros(1,total_n)
for i=1:size(opt_value.Pg)[2]
    time_stamp[i]=i
end

Solar=[time_stamp; opt_value.Pg; opt_value.Qg]';
# println(size(Solar))
CSV.write("solar.csv", DataFrame(Solar, [:timestamp, :Pg_1, :Pg_2, :Pg_3,
    :Pg_4, :Pg_5, :Pg_6, :Pg_7, :Pg_8, :Pg_9, :Pg_10, :Pg_11, :Pg_12, :Qg_1,
    :Qg_2, :Qg_3, :Qg_4, :Qg_5, :Qg_6, :Qg_7, :Qg_8, :Qg_9, :Qg_10,
    :Qg_11, :Qg_12]));

# write the battery file
Battery=[time_stamp; opt_value.R; opt_value.C; opt_value.B]';
# println(size(Battery))
CSV.write("battery.csv", DataFrame(Battery, [:timestamp, :R_1, :R_2, :R_3,
    :R_4, :R_5, :R_6, :R_7, :R_8, :R_9, :R_10, :R_11, :R_12, :C_1,
    :C_2, :C_3, :C_4, :C_5, :C_6, :C_7, :C_8, :C_9, :C_10,
    :C_11, :C_12, :B_1, :B_2, :B_3, :B_4, :B_5, :B_6, :B_7, :B_8, :B_9, :B_10,
     :B_11, :B_12]));

println("GML - Finish writting files! ")
end
