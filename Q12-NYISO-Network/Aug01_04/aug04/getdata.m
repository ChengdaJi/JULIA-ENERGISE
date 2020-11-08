initial_pos = 1;
end_pos = 288;
opt_length = end_pos-initial_pos+1;

% 1 for bat only
total_cost_1 = 0;
p0_1 = zeros(opt_length,1);
pg_1 = zeros(opt_length,1);
R_1 = zeros(opt_length,1);
prsrv_1 = zeros(opt_length,1);
B_1  = zeros(opt_length,1);
pcul_1 = zeros(opt_length,1);
st_1  = zeros(opt_length,1);
pd_1 = zeros(opt_length,1);

cd('GML_Bat') 
cd('result') 
for i = 1:opt_length
    s1 = 'M1P_rate';
%     s2 = ratelist(j);
    s3 = 'Time';
    s4 = string(initial_pos+i-1);
    s5 = '.csv';

    path = strcat(s3,s4,s5);
    T = readtable(path);

    temp_cost = table2array(T(1,1));
    total_cost_1(i) = temp_cost;
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
    elseif temp_statue == "SLOW_PROGRESS"
        i
    end
end
cd('..') 
cd('..') 


% 2 for  peak shaving
total_cost_2 = 0;
p0_2 = zeros(opt_length,1);
pg_2 = zeros(opt_length,1);
R_2 = zeros(opt_length,1);
prsrv_2 = zeros(opt_length,1);
B_2  = zeros(opt_length,1);
pcul_2 = zeros(opt_length,1);
st_2  = zeros(opt_length,1);

cd('GML_Bat_Peak') 
cd('result') 
for i = 1:opt_length
    s1 = 'M1P_rate';
%     s2 = ratelist(j);
    s3 = 'Time';
    s4 = string(initial_pos+i-1);
    s5 = '.csv';

    path = strcat(s3,s4,s5);
    T = readtable(path);

    temp_cost = table2array(T(1,1));
    total_cost_2(i) = temp_cost;
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
    elseif temp_statue == "SLOW_PROGRESS"
        i
    end
end
cd('..') 
cd('..') 


% 3 for  peak shaving
total_cost_3 = 0;
p0_3 = zeros(opt_length,1);
pg_3 = zeros(opt_length,1);
R_3 = zeros(opt_length,1);
prsrv_3 = zeros(opt_length,1);
B_3  = zeros(opt_length,1);
pcul_3 = zeros(opt_length,1);
st_3  = zeros(opt_length,1);
cd('GML_No') 
cd('result') 
for i = 1:opt_length
    s1 = 'M1P_rate';
%     s2 = ratelist(j);
    s3 = 'Time';
    s4 = string(initial_pos+i-1);
    s5 = '.csv';

    path = strcat(s3,s4,s5);
    T = readtable(path);

    temp_cost = table2array(T(1,1));
    total_cost_3(i) = temp_cost;
%         Cost_real_time = Cost_real_time + table2array(T(1,11));
%         Cost_solar_cul = Cost_solar_cul + table2array(T(1,12));
%         Revenue_rsrv = Revenue_rsrv + table2array(T(1,13));


    temp_pg = table2array(T(1,3));
    pg_3(i) = temp_pg;

    temp_R = table2array(T(1,5));
    R_3(i) = temp_R;

    temp_p0 = table2array(T(1,7));
    p0_3(i) = temp_p0;
    
    temp_pd = table2array(T(1,6));
    pd_3(i) = temp_pd;

    temp_prsrv = table2array(T(1,8));
    prsrv_3(i) = temp_prsrv;

    temp_b = table2array(T(1,4));
    B_3(i) = temp_b;

    temp_st = table2array(T(1,2));
    st_3(i) = temp_st;
    
    temp_pt = table2array(T(1,9));
    pcul_3(i) = temp_pt;
    
    temp_statue = table2array(T(1,10));

    if temp_statue == "INFEASIBLE"
        i
    elseif temp_statue == "SLOW_PROGRESS"
        i
    end
end
cd('..') 
cd('..') 

s = 'trace_aug04_2019.csv';
path = strcat(s);
T = readtable(path);
lambda_da = table2array(T(:,5));
lambda_rt = table2array(T(:,4));
lambda_rt = lambda_rt(1:end);
lambda_da = lambda_da(1:end);

rt_cost_bat_GML = total_cost_1;
rt_cost_bat_peak_GML = total_cost_2;
rt_cost_no_GML = total_cost_3;

rt_cost_bat = vpa(lambda_rt(initial_pos+1:end_pos+1)' * p0_1(initial_pos:end_pos) *100/12, 7)
rt_cost_bat_peak = vpa(lambda_rt(initial_pos+1:end_pos+1)' * p0_2(initial_pos:end_pos) *100/12, 7)
rt_cost_no = vpa(lambda_rt(initial_pos+1:end_pos+1)' * p0_3(initial_pos:end_pos) *100/12, 7)


peak_bat = vpa(max(p0_1(initial_pos:end_pos))*1000000, 7)
peak_bat_peak = vpa(max(p0_2(initial_pos:end_pos))*1000000, 7)
peak_no = vpa(max(p0_3(initial_pos:end_pos))*1000000, 7)


sol_cul_bat = vpa(sum(pcul_1(initial_pos:end_pos))*15*100/12, 7)
sol_cul_bat_peak = vpa(sum(pcul_2(initial_pos:end_pos))*15*100/12, 7)
sol_cul_no = vpa(sum(pcul_3(initial_pos:end_pos))*15*100/12, 7)

