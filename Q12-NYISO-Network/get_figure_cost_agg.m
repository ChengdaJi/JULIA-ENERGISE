s = 'SEP17 noise/trace_aug01_2019.csv';
path = strcat(s);
T = readtable(path);
lambda_da = table2array(T(:,5));
lambda_rt = table2array(T(:,4));
lambda_rt = lambda_rt(1:288);
lambda_da = lambda_da(1:288);
% 1 for no peak shaving
cost_agg_1 = zeros(288,1);


cd('SEP17_CUL_NOISE/')
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

    temp_cost_agg = table2array(T(1,11));
    cost_agg_1(i) = temp_cost_agg;

end
cd('..') 
cd('..') 

% 2 for with peak shaving
cost_agg_2 = zeros(288,1);


% cd('SEP17_CUL_NOISE/')
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

    temp_cost_agg = table2array(T(1,11));
    cost_agg_2(i) = temp_cost_agg;

end
cd('..') 
cd('..')
cd('..')

% 2 for with peak shaving
cost_agg_3 = zeros(288,1);

cd('SEP17_CUL_NOISE_Bat/')
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

    temp_cost_agg = table2array(T(1,11));
    cost_agg_3(i) = temp_cost_agg;

end

times = 1:288;
figure('Name',"p0`")
hold on
% yyaxis right
% l = plot(times,lambda_rt,'k','linestyle','--','linewidth',3);
% set(gca,'Ycolor','black','FontSize',20)
% ylabel('[$\$/$MWHr]','Interpreter','latex','FontSize',20,'Color','black');
% ylim([0,50])


% plot(times,cost_agg_1,'r','linestyle','-','linewidth',3);
% plot(times,cost_agg_2,'b','linestyle','-','linewidth',3);
% plot(times,cost_agg_3,'k','linestyle','-.','linewidth',3);

% bar(cost_agg_3-2*10^5);
% 
% bar(cost_agg_2-2*10^5);
% 
% bar(cost_agg_1-2*10^5);
yyaxis left
bar(cost_agg_3-cost_agg_2, );

bar(cost_agg_2-cost_agg_1);


% bar(cost_agg_1-2*10^5);

set(gca,'Ycolor','black','FontSize',20)
legend('No shaving and battery - With battery and shaving','With battery and shaving - With battery only')

ylabel('Cost difference','Interpreter','latex','FontSize',20,'Color','black');
set(gca,'xtick',[1,36,72,108,144,180,216,252,288], 'xticklabel',{'12AM','3AM','6AM','9AM','12PM','3PM','6PM','9PM','12AM'});

yyaxis right
l = plot(times,lambda_rt,'k','linestyle','--','linewidth',3);
set(gca,'Ycolor','black','FontSize',20)
ylabel('[$\$/$MWHr]','Interpreter','latex','FontSize',20,'Color','black');
ylim([0,50])
xlim([1,288]);
% ylim([-20,40]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 35, 6], 'PaperUnits', 'Inches', 'PaperSize', [35, 6])
% plot([231,231],[-20,40],'linestyle','--','linewidth',3)
grid on;
saveas(gcf,'P','epsc')
cd('../../..')