B_cap=[0 10 20 30];
i=1;
for B = 0:10:30
    
    filename = strcat('B=',string(B),'MWh.csv');
    T = readtable(filename);
    p_gen{i} = table2array(T(1:4,1));
    q_gen{i} = table2array(T(1:4,2));
    v{i} = table2array(T(1:4,3));
    theta{i} = table2array(T(1:4,4));
    p_gen_star{i} = table2array(T(1:4,5));
    q_gen_star{i} = table2array(T(1:4,6));
    v_star{i} = table2array(T(1:4,7));
    theta_star{i} = table2array(T(1:4,8));
    
    i=i+1;
end
    
% voltage difference
figure(1)
hold on
plot(1:4, v{1}-v_star{1},'b','linestyle','-','linewidth',3);
% plot(1:4, v_star{1},'b','linestyle',':','linewidth',3);
plot(1:4, v{2}-v_star{2},'r','linestyle','-','linewidth',3);
% plot(1:4, v_star{2},'r','linestyle',':','linewidth',3);
plot(1:4, v{3}-v_star{3},'g','linestyle','-','linewidth',3);
% plot(1:4, v_star{3},'g','linestyle',':','linewidth',3);
plot(1:4, v{4}-v_star{4},'y','linestyle','-','linewidth',3);
% plot(1:4, v_star{4},'y','linestyle',':','linewidth',3);
title('voltage')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('0','10','20','30')

figure(2)
hold on
plot(1:4, theta{1}-theta_star{1},'b','linestyle','-','linewidth',3);
% plot(1:4, theta_star{1},'b','linestyle',':','linewidth',3);
plot(1:4, theta{2}-theta_star{2},'r','linestyle','-','linewidth',3);
% plot(1:4, theta_star{2},'r','linestyle',':','linewidth',3);
plot(1:4, theta{3}-theta_star{3},'g','linestyle','-','linewidth',3);
% plot(1:4, theta_star{3},'g','linestyle',':','linewidth',3);
plot(1:4, theta{4}-theta_star{4},'y','linestyle','-','linewidth',3);
% plot(1:4, v_star{4},'y','linestyle',':','linewidth',3);
title('theta')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('0','10','20','30')

figure(3)
hold on
plot(1:4, p_gen{1}-p_gen_star{1},'b','linestyle','-','linewidth',3);
% plot(1:4, p_gen_star{1},'b','linestyle',':','linewidth',3);
plot(1:4, p_gen{2}-p_gen_star{2},'r','linestyle','-','linewidth',3);
% plot(1:4, p_gen_star{2},'r','linestyle',':','linewidth',3);
plot(1:4, p_gen{3}-p_gen_star{3},'g','linestyle','-','linewidth',3);
% plot(1:4, p_gen_star{3},'g','linestyle',':','linewidth',3);
plot(1:4, p_gen{4}-p_gen_star{4},'y','linestyle','-','linewidth',3);
% plot(1:4, p_gen_star{4},'y','linestyle',':','linewidth',3);
title('p_gen')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('0','10','20','30')

figure(4)
hold on
plot(1:4, q_gen{1}-q_gen_star{1},'b','linestyle','-','linewidth',3);
% plot(1:4, q_gen_star{1},'b','linestyle',':','linewidth',3);
plot(1:4, q_gen{2}-q_gen_star{2},'r','linestyle','-','linewidth',3);
% plot(1:4, q_gen_star{2},'r','linestyle',':','linewidth',3);
plot(1:4, q_gen{3}-q_gen_star{3},'g','linestyle','-','linewidth',3);
% plot(1:4, q_gen_star{3},'g','linestyle',':','linewidth',3);
plot(1:4, q_gen{4}-q_gen_star{4},'y','linestyle','-','linewidth',3);
% plot(1:4, q_gen_star{4},'y','linestyle',':','linewidth',3);
title('q_gen')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('0','10','20','30')


% figure(2)
% hold on
% plot(1:4, theta,'b','linestyle','--','linewidth',3);
% plot(1:4, theta_star,'r','linestyle',':','linewidth',3);
% title('THETA')
% xlabel('bus')
% ylabel('p.u.')
% set(gca,'Ycolor','black','FontSize',20)
% legend('GML result','PowerModel')
% 
% figure(3)
% hold on
% plot(1:4, P_gen,'b','linestyle','--','linewidth',3);
% plot(1:4, P_gen_star,'r','linestyle',':','linewidth',3);
% title('Real Gen')
% xlabel('bus')
% ylabel('p.u.')
% set(gca,'Ycolor','black','FontSize',20)
% legend('GML result','PowerModel')
% 
% figure(4)
% hold on
% plot(1:4, Q_gen,'b','linestyle','--','linewidth',3);
% plot(1:4, Q_gen_star,'r','linestyle',':','linewidth',3);
% title('Reactive Gen')
% xlabel('bus')
% ylabel('p.u.')
% set(gca,'Ycolor','black','FontSize',20)
% legend('GML result','PowerModel')