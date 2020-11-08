% s = 'P_gen.xlsx';
% path = strcat(s);
% T = readtable(path);
% pd = table2array(T(:,2));
% pd = p0_3
% lambda_rt = table2array(T(:,4));
% lambda_rt = lambda_rt(1:288);
% lambda_da = lambda_da(1:288);


times = 1:288;
figure('Name',"p0")
hold on
yyaxis right
l = plot(times,lambda_rt(2:289),'k','linestyle','--','linewidth',3);
set(gca,'Ycolor','black','FontSize',20)
ylabel('[$\$/$MWHr]','Interpreter','latex','FontSize',20,'Color','black');
ylim([0,50])

yyaxis left
plot(times,(p0_1)*100,'r','linestyle','-','linewidth',3);
plot(times,(p0_2)*100,'b','linestyle',':','linewidth',3);
plot(times,(p0_3)*100,'k','linestyle','-.','linewidth',3);


set(gca,'Ycolor','black','FontSize',20)
legend('$P_{0}$ no shaving','$P_{0}$ with shaving','$P_{d}$','LMP','Interpreter','latex')

ylabel('$P$\,\,[MW]','Interpreter','latex','FontSize',20,'Color','black');
set(gca,'xtick',[1,36,72,108,144,180,216,252,288], 'xticklabel',{'12AM','3AM','6AM','9AM','12PM','3PM','6PM','9PM','12AM'});
xlim([1,288]);
% ylim([-20,40]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 35, 6], 'PaperUnits', 'Inches', 'PaperSize', [35, 6])
% plot([231,231],[-20,40],'linestyle','--','linewidth',3)
grid on;
saveas(gcf,'P','epsc')


figure('Name',"R")
hold on
yyaxis right
l = plot(times,lambda_rt(2:289),'k','linestyle','--','linewidth',3);
set(gca,'Ycolor','black','FontSize',20)
ylabel('[$\$/$MWHr]','Interpreter','latex','FontSize',20,'Color','black');
ylim([0,50])

yyaxis left
plot(times,(R_1)*100,'r','linestyle','-','linewidth',3);
plot(times,(R_2)*100,'b','linestyle','-','linewidth',3);
ylim([-75,75])

set(gca,'Ycolor','black','FontSize',20)
legend('$R$ no shaving','$R$ with shaving','LMP','Interpreter','latex')

ylabel('$R$\,\,[MW]','Interpreter','latex','FontSize',20,'Color','black');
set(gca,'xtick',[1,36,72,108,144,180,216,252,288], 'xticklabel',{'12AM','3AM','6AM','9AM','12PM','3PM','6PM','9PM','12AM'});
xlim([1,288]);
% ylim([-20,40]);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 35, 6], 'PaperUnits', 'Inches', 'PaperSize', [35, 6])
saveas(gcf,'R','epsc')
% plot([231,231],[-20,40],'linestyle','--','linewidth',3)
grid on;