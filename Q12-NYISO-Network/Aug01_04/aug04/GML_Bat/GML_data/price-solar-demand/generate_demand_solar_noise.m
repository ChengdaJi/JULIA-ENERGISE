% demand = load('normalized_demand.mat','normalized_demand');
% pd=demand.normalized_demand;
% sum(pd,1)
% figure(1)
% plot(1:576,pd(:,1))
% % hold on 
% plot(1:576,pd(:,2));

% for day = 1:4
%     for time = 1:288
%         
%         outputname = strcat('demand_noise/demand_noise_aug0',...
%              num2str(day),'_time',num2str(time),'.mat')
%         demand_noise = zeros(240,24);
%         predict_time = 1;
%         for predict_std = 0.01:0.005/23:0.015
%             demand_noise(:,predict_time)=...
%                 randn(240,1)*predict_std;
%             predict_time = predict_time+1;
%         end
%         save(outputname, 'demand_noise')
%     end
% end


for day = 1:4
    for time = 1:288
        
        outputname = strcat('solar_noise/solar_noise_aug0',...
             num2str(day),'_time',num2str(time),'.mat')
        solar_noise = zeros(240,24);
        predict_time = 1;
        for predict_std = 0.1:0.05/23:0.15
            solar_noise(:, predict_time)=...
                randn(240,1)*predict_std;
            predict_time = predict_time+1;
        end
        save(outputname, 'solar_noise')
    end
end