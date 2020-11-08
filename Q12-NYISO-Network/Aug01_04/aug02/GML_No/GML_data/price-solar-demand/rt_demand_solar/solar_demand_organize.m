%% organize day_ahead demand
% for i=1:4
%     name1 = strcat('aug0',num2str(i),'_da.mat');
%     demand1 = load(name);
%     name2 = strcat('aug0',num2str(i),'_da.mat')
%     demand2 = load(name);
%     day_ahead = [demand1.aug1; demand1.aug1]/2;
%     day_ahead = demand1.aug1;
%     outputname = strcat('demand_da_aug0',num2str(i),'.mat')
%     save(outputname, 'day_ahead')
% end


% da = load('demand_da_aug01.mat','day_ahead');
% sum(da.day_ahead, 1)
%% organize real time demand
% for i=1:4
%     name = strcat('aug0',num2str(i),'_realtime.mat')
%     demand = load(name);
%     real_time = demand.aug1;
%     outputname = strcat('demand_aug0',num2str(i),'.mat')
%     save(outputname, 'real_time')
% end

%% organize solar real time demand
% for i=1:4
%     name = strcat('aug0',num2str(i),'_solar.mat')
%     solar = load(name);
%     real_time = solar.result;
%     outputname = strcat('solar_aug0',num2str(i),'.mat')
%     save(outputname, 'real_time')
% end

%% organize solar real time demand
for i=1:4
    name = strcat('aug0',num2str(i),'_solar.mat')
    solar = load(name);
    real_time= solar.result;
    day_ahead = zeros(576,240);
    for t=1:576
        group = ((ceil(t/12)-1)*12+1);
        day_ahead (t,:) = real_time(group,:);
    end
    outputname = strcat('solar_aug0',num2str(i),'_day_ahead.mat')
    save(outputname, 'day_ahead')
end
