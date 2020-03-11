% demand_data_raw=readmatrix("demand_feb_6_7.csv");
% original_nonunified = demand_data_raw(1:576,3:14);
% for or_n=1:12
%     for number=1:20
%         temp_rand = 0.05*randn(576,1);
%         normalized_demand(:, (or_n-1)*20+number) = abs((original_nonunified(:, or_n)+temp_rand))/sum(abs(original_nonunified(:, or_n)+temp_rand));
%     end
% end
% save('normalized_demand.mat','normalized_demand')

% data=load('normalized_demand.mat','normalized_demand');
% demand=data.normalized_demand;
% normalized_demand_da=zeros(576,240);
% for i=1:240
%     for t=1:576
%         normalized_demand_da(t,i)=demand((ceil(t/12)-1)*12+1,i);
%     end
% end
% figure (1)
% hold on
% 
% for i=1:120
%     plot(1:576, normalized_demand_da(:,i))
% end
% 
% save('normalized_demand_da.mat','normalized_demand_da')


data=load('normalized_solar.mat','normalized_solar');
solar=data.normalized_solar;
normalized_solar_da=zeros(576,240);
for i=1:240
    abs(1+0.1*randn(1))
    for t=1:576
        if mod(t,12)==1
            mult = abs(1+0.1*randn(1));
        end
        normalized_solar_da(t,i)=mult*solar((ceil(t/12)-1)*12+1,i);
%         normalized_solar_da2(t,i)=solar((ceil(t/12)-1)*12+1,i);
    end
end
figure (1)
hold on

for i=1:120
    plot(1:576, normalized_solar_da(:,i))
end

save('normalized_solar_da.mat','normalized_solar_da')



% demand_solar_raw=readmatrix("solar.csv");
% original_solar = demand_solar_raw(:,2:13);
% normalized_solar=zeros(576,240);
% for or_n=1:12
%     for number=1:20
%         temp_rand = 0.25*randn(576,1);
%         normalized_solar(:, (or_n-1)*20+number) = ...
%             abs((original_solar(:, or_n)+original_solar(:, or_n).*temp_rand))/...
%             sum(abs(original_solar(:, or_n)+original_solar(:, or_n).*temp_rand));
%     end
% end
% figure (1)
% hold on
% for i=1:240
% plot(1:576, normalized_solar(:,i))
% end
% 
% save('normalized_solar.mat','normalized_solar')