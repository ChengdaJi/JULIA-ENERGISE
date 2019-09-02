% This file is created to process the csv file from NYISO
% This file generate MAT file for matlab
% Chengda Ji
% Jun 3, 2019

clear all;
clc;
% addpath('data_raw');
selected_zone='HUD VL';

%% LMP real-time data
lambda_rt = [];
T_LMP_RT = readtable('data_raw/energy/20190529realtime_zone.csv');
T_LMP_RT.Properties.VariableNames={'Time','Name','PTID','LMP','MCL','MCC'};

for T_LMP_RT_rownumber = 1:length(T_LMP_RT.Name)
    if strcmp(T_LMP_RT.Name{T_LMP_RT_rownumber},selected_zone)
        T_LMP_RT.LMP(T_LMP_RT_rownumber);
        lambda_rt=[lambda_rt,T_LMP_RT.LMP(T_LMP_RT_rownumber)];
    end
end
T_LMP_RT = readtable('data_raw/energy/20190530realtime_zone.csv');
T_LMP_RT.Properties.VariableNames={'Time','Name','PTID','LMP','MCL','MCC'};

for T_LMP_RT_rownumber = 1:length(T_LMP_RT.Name)
    if strcmp(T_LMP_RT.Name{T_LMP_RT_rownumber},selected_zone)
        T_LMP_RT.LMP(T_LMP_RT_rownumber);
        lambda_rt=[lambda_rt,T_LMP_RT.LMP(T_LMP_RT_rownumber)];
    end
end
size(lambda_rt)

%% LMP Day-ahead  
lambda_da_raw = [];
T_LMP_DA = readtable('data_raw/energy/20190529damlbmp_zone.csv');
T_LMP_DA.Properties.VariableNames={'Time','Name','PTID','LMP','MCL','MCC'};

for T_LMP_RT_rownumber = 1:length(T_LMP_DA.Name)
    if strcmp(T_LMP_DA.Name{T_LMP_RT_rownumber},selected_zone)
        T_LMP_DA.LMP(T_LMP_RT_rownumber);
        lambda_da_raw=[lambda_da_raw,T_LMP_DA.LMP(T_LMP_RT_rownumber)];
    end
end
T_LMP_DA = readtable('data_raw/energy/20190530damlbmp_zone.csv');
T_LMP_DA.Properties.VariableNames={'Time','Name','PTID','LMP','MCL','MCC'};

for T_LMP_RT_rownumber = 1:length(T_LMP_DA.Name)
    if strcmp(T_LMP_DA.Name{T_LMP_RT_rownumber},selected_zone)
        T_LMP_DA.LMP(T_LMP_RT_rownumber);
        lambda_da_raw=[lambda_da_raw,T_LMP_DA.LMP(T_LMP_RT_rownumber)];
    end
end
T_LMP_DA = readtable('data_raw/energy/20190531damlbmp_zone.csv');
T_LMP_DA.Properties.VariableNames={'Time','Name','PTID','LMP','MCL','MCC'};

for T_LMP_RT_rownumber = 1:length(T_LMP_DA.Name)
    if strcmp(T_LMP_DA.Name{T_LMP_RT_rownumber},selected_zone)
        T_LMP_DA.LMP(T_LMP_RT_rownumber);
        lambda_da_raw=[lambda_da_raw,T_LMP_DA.LMP(T_LMP_RT_rownumber)];
        break
    end
end
lambda_da=[];
for lambda_da_number=1:length(lambda_da_raw)-1
    lambda_da_inter = interp1([lambda_da_number,lambda_da_number+1],...
        [lambda_da_raw(lambda_da_number),...
        lambda_da_raw(lambda_da_number+1)],...
        lambda_da_number:1/12: lambda_da_number+1-1/12);
    lambda_da = [lambda_da, lambda_da_inter];
end





%% Ancillary real-time data
alpha_rt_10 = [];
alpha_rt_30 = [];
T_ANL_RT = readtable('data_raw/ancillary/20190529rtasp.csv');
T_ANL_RT.Properties.VariableNames={'Time','TimeZone','Name','PTID',...
    'Spin10','Nonsyn10','Opr30','RegCap','RegMov'};

for T_ANL_RL_rownumber = 1:length(T_ANL_RT.Name)
    if strcmp(T_ANL_RT.Name{T_ANL_RL_rownumber},selected_zone)
        alpha_rt_10=[alpha_rt_10,T_ANL_RT.Spin10(T_ANL_RL_rownumber)];
        alpha_rt_30=[alpha_rt_30,T_ANL_RT.Opr30(T_ANL_RL_rownumber)];
    end
end
T_ANL_RT = readtable('data_raw/ancillary/20190530rtasp.csv');
T_ANL_RT.Properties.VariableNames={'Time','TimeZone','Name','PTID',...
    'Spin10','Nonsyn10','Opr30','RegCap','RegMov'};

for T_ANL_RL_rownumber = 1:length(T_ANL_RT.Name)
    if strcmp(T_ANL_RT.Name{T_ANL_RL_rownumber},selected_zone)
        alpha_rt_10=[alpha_rt_10,T_ANL_RT.Spin10(T_ANL_RL_rownumber)];
        alpha_rt_30=[alpha_rt_30,T_ANL_RT.Opr30(T_ANL_RL_rownumber)];
    end
end
size(alpha_rt_10)
size(alpha_rt_30)



%% Ancillary Day-ahead  
alpha_da_raw_10 = [];
alpha_da_raw_30 = [];
T_ANL_DA = readtable('data_raw/ancillary/20190529damasp.csv');
T_ANL_DA.Properties.VariableNames={'Time','TimeZone','Name','PTID',...
    'Spin10','Nonsyn10','Opr30','RegCap'};

for T_ANL_DA_rownumber = 1:length(T_ANL_DA.Name)
    if strcmp(T_ANL_DA.Name{T_ANL_DA_rownumber},selected_zone)
        alpha_da_raw_10=[alpha_da_raw_10,T_ANL_DA.Spin10(T_ANL_DA_rownumber)];
        alpha_da_raw_30=[alpha_da_raw_30,T_ANL_DA.Opr30(T_ANL_DA_rownumber)];
    end
end
T_ANL_DA = readtable('data_raw/ancillary/20190530damasp.csv');
T_ANL_DA.Properties.VariableNames={'Time','TimeZone','Name','PTID',...
    'Spin10','Nonsyn10','Opr30','RegCap'};

for T_ANL_DA_rownumber = 1:length(T_ANL_DA.Name)
    if strcmp(T_ANL_DA.Name{T_ANL_DA_rownumber},selected_zone)
        alpha_da_raw_10=[alpha_da_raw_10,T_ANL_DA.Spin10(T_ANL_DA_rownumber)];
        alpha_da_raw_30=[alpha_da_raw_30,T_ANL_DA.Opr30(T_ANL_DA_rownumber)];
    end
end
T_ANL_DA = readtable('data_raw/ancillary/20190531damasp.csv');
T_ANL_DA.Properties.VariableNames={'Time','TimeZone','Name','PTID',...
    'Spin10','Nonsyn10','Opr30','RegCap'};

for T_ANL_DA_rownumber = 1:length(T_ANL_DA.Name)
    if strcmp(T_ANL_DA.Name{T_ANL_DA_rownumber},selected_zone)
        alpha_da_raw_10=[alpha_da_raw_10,T_ANL_DA.Spin10(T_ANL_DA_rownumber)];
        alpha_da_raw_30=[alpha_da_raw_30,T_ANL_DA.Opr30(T_ANL_DA_rownumber)];
        break
    end
end
size(alpha_da_raw_10)
size(alpha_da_raw_30)


alpha_da_10=[];
alpha_da_30=[];
for alpha_da_number=1:length(alpha_da_raw_10)-1
    alpha_da_inter_10 = interp1([alpha_da_number,alpha_da_number+1],...
        [alpha_da_raw_10(alpha_da_number),...
        alpha_da_raw_10(alpha_da_number+1)],...
        alpha_da_number:1/12: alpha_da_number+1-1/12);
    alpha_da_inter_30 = interp1([alpha_da_number,alpha_da_number+1],...
        [alpha_da_raw_30(alpha_da_number),...
        alpha_da_raw_30(alpha_da_number+1)],...
        alpha_da_number:1/12: alpha_da_number+1-1/12);
    alpha_da_10 = [alpha_da_10, alpha_da_inter_10];
    alpha_da_30 = [alpha_da_30, alpha_da_inter_30];
end
size(alpha_da_10)
size(alpha_da_30)
alpha_da_10



save('price_data.mat','lambda_rt','lambda_da','alpha_rt_10',...
    'alpha_rt_30','alpha_da_10','alpha_da_30');





