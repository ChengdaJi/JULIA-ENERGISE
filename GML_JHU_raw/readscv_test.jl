using CSV
using DataFrames
using MAT

# using JuliaDB
# Bank1_raw=CSV.read("demand_and_solar_allendale_bank139_aug2016.csv");
# Bank2_raw=CSV.read("demand_and_solar_allendale_bank239_aug2016.csv");
# Bank3_raw=CSV.read("demand_and_solar_blooming_grove_bank276_aug2016.csv");
# Bank1_raw_5th=zeros(288,4)
# Bank2_raw_5th=zeros(288,4)
# Bank3_raw_5th=zeros(288,4)
# for i=1:288
#     k=(i-1)*5+1;
#     Bank1_raw_5th[i,:]=[Bank1_raw[k,4],Bank1_raw[k,9],Bank1_raw[k,14],
#     Bank1_raw[k,19]]
#     Bank2_raw_5th[i,:]=[Bank2_raw[k,4],Bank2_raw[k,9],Bank2_raw[k,14],
#     Bank2_raw[k,19]]
#     Bank3_raw_5th[i,:]=[Bank3_raw[k,4],Bank3_raw[k,9],Bank3_raw[k,14],
#     Bank3_raw[k,19]]
# end
# Pd=hcat(Bank1_raw_5th/1000, Bank2_raw_5th/1000, Bank3_raw_5th/1000)
# CSV.write("Pd_Feb1.csv",DataFrame(Pd))
# Pd_2=CSV.read("Pd_Feb1.csv");

# Pg1_raw=CSV.read("02.2017_data_al.csv");
# Pg2_raw=CSV.read("02.2017_data_bl.csv");
#
# Pg_trend=zeros(288,2)
#
# for i=1:288
#     k=(i-1)*5+2;
#     Pg_trend[i,:]=[Pg1_raw[k,8],Pg2_raw[k,8]]
# end
#
# println(size(Pg_trend))
# CSV.write("Pg_trend_Feb1.csv",DataFrame(Pg_trend))
Price_mat=matread("data/price.mat");
Alpha_raw=Price_mat["alpha_10"][1:288];
Lambda_raw=Price_mat["Lambda"][1:288];
Lambda=zeros(288,1)
Alpha=zeros(288,1)
Lambda[:,1]=Lambda_raw
Alpha[:,1]=Alpha_raw
println(size(Lambda))
CSV.write("Alpha_Feb1.csv",DataFrame(Alpha))
CSV.write("Lambda_Feb1.csv",DataFrame(Lambda))
