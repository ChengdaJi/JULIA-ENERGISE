% demand_noise = {};
% for feeder = 1:12
%     normalized_noise = randn(288,24);
%     demand_noise_one = sqrt(0.01:0.01/23:0.02);
%     demand_noise_var = repmat(demand_noise_one, 288,1);
%     demand_noise{feeder} = normalized_noise.*demand_noise_var;
% end
% save('demand_noise.mat','demand_noise')

% solar_noise = {};
% dest = 0.1;
% for feeder = 1:12
%     normalized_noise = randn(288,24);
%     solar_noise_one = sqrt(0.01:(dest-0.01)/23:dest);
%     solar_noise_var = repmat(solar_noise_one, 288,1);
%     solar_noise{feeder} = normalized_noise.*solar_noise_var;
% end
% save('solar_noise_01.mat','solar_noise')

pd = load('solar_noise_0025.mat');
var(pd.solar_noise{1}(:,1))
var(pd.solar_noise{1}(:,24))