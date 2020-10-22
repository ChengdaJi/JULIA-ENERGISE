% 1 for no peak shaving
total_cost_1 = 0;
p0_1 = zeros(288,1);
pg_1 = zeros(288,1);
R_1 = zeros(288,1);
prsrv_1 = zeros(288,1);
B_1  = zeros(288,1);
pcul_1 = zeros(288,1);
st_1  = zeros(288,1);
pd_1 = zeros(288,1);

cd('GML_JHU_NYISO_NO_peak') 
cd('result') 
for i = 1:288
    s1 = 'M1P_rate';
%     s2 = ratelist(j);
    s3 = 'Time';
    s4 = string(i);
    s5 = '.csv';

    path = strcat(s3,s4,s5);
    T = readtable(path);

    temp_cost = table2array(T(1,1));
    total_cost_1 = total_cost_1 + temp_cost;
%         Cost_real_time = Cost_real_time + table2array(T(1,11));
%         Cost_solar_cul = Cost_solar_cul + table2array(T(1,12));
%         Revenue_rsrv = Revenue_rsrv + table2array(T(1,13));


    temp_pg = table2array(T(1,3));
    pg_1(i) = temp_pg;

    temp_R = table2array(T(1,5));
    R_1(i) = temp_R;

    temp_p0 = table2array(T(1,7));
    p0_1(i) = temp_p0;
    
    temp_pd = table2array(T(1,6));
    pd_1(i) = temp_pd;

    temp_prsrv = table2array(T(1,8));
    prsrv_1(i) = temp_prsrv;

    temp_b = table2array(T(1,4));
    B_1(i) = temp_b;

    temp_st = table2array(T(1,2));
    st_1(i) = temp_st;
    
    temp_pt = table2array(T(1,9));
    pcul_1(i) = temp_pt;
    
    temp_statue = table2array(T(1,10));
%     print("no peak shaving")
    if temp_statue == "INFEASIBLE"
        i
    end
end
cd('..') 
cd('..') 


% 2 for  peak shaving
total_cost_2 = 0;
p0_2 = zeros(288,1);
pg_2 = zeros(288,1);
R_2 = zeros(288,1);
prsrv_2 = zeros(288,1);
B_2  = zeros(288,1);
pcul_2 = zeros(288,1);
st_2  = zeros(288,1);

cd('GML_JHU_NYISO_W_peak') 
cd('result') 
for i = 1:288
    s1 = 'M1P_rate';
%     s2 = ratelist(j);
    s3 = 'Time';
    s4 = string(i);
    s5 = '.csv';

    path = strcat(s3,s4,s5);
    T = readtable(path);

    temp_cost = table2array(T(1,1));
    total_cost_2 = total_cost_2 + temp_cost;
%         Cost_real_time = Cost_real_time + table2array(T(1,11));
%         Cost_solar_cul = Cost_solar_cul + table2array(T(1,12));
%         Revenue_rsrv = Revenue_rsrv + table2array(T(1,13));


    temp_pg = table2array(T(1,3));
    pg_2(i) = temp_pg;

    temp_R = table2array(T(1,5));
    R_2(i) = temp_R;

    temp_p0 = table2array(T(1,7));
    p0_2(i) = temp_p0;
    
    temp_pd = table2array(T(1,6));
    pd_2(i) = temp_pd;

    temp_prsrv = table2array(T(1,8));
    prsrv_2(i) = temp_prsrv;

    temp_b = table2array(T(1,4));
    B_2(i) = temp_b;

    temp_st = table2array(T(1,2));
    st_2(i) = temp_st;
    
    temp_pt = table2array(T(1,9));
    pcul_2(i) = temp_pt;
    
    temp_statue = table2array(T(1,10));

    if temp_statue == "INFEASIBLE"
        i
    end
end
cd('..') 
cd('..') 

s = 'trace_aug01_2019.csv';
path = strcat(s);
T = readtable(path);
lambda_da = table2array(T(:,5));
lambda_rt = table2array(T(:,4));
lambda_rt = lambda_rt(1:288);
lambda_da = lambda_da(1:288);

initial_slot =17*12+1;
end_slot = 20*12;

peak_no = max(p0_1(initial_slot:end_slot))*1000000
peak_with = max(p0_2(initial_slot:end_slot))*1000000
% peak_peak = max(p0_3(initial_slot:end_slot))*1000000


rt_cost_no = lambda_rt(initial_slot:end_slot)' * p0_1(initial_slot:end_slot) *100/12
rt_cost_with = lambda_rt(initial_slot:end_slot)' * p0_2(initial_slot:end_slot) *100/12
% rt_cost_peak = lambda_rt(initial_slot:end_slot)' * p0_3(initial_slot:end_slot) *100/12


sol_cul_1 = sum(pcul_1(initial_slot:end_slot))*15*100/12
sol_cul_2 = sum(pcul_2(initial_slot:end_slot))*15*100/12
