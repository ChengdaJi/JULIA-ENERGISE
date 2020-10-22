% B_cap=[0 10 20 30];
i=1;
for Q = 0:0.25:1
    
    filename = strcat('Q=',string(Q),'_P_gen.csv');
    T = readtable(filename);
    p_gen{i} = table2array(T(:,1));
    p_gen_star{i} = table2array(T(:,2));
    
    filename = strcat('Q=',string(Q),'_Q_gen.csv');
    T = readtable(filename);
    q_gen{i} = table2array(T(:,1));
    q_gen_star{i} = table2array(T(:,2));
    filename = strcat('Q=',string(Q),'_v.csv');
    T = readtable(filename);
    v{i} = table2array(T(:,1));
    v_star{i} = table2array(T(:,2));
    filename = strcat('Q=',string(Q),'_theta.csv');
    T = readtable(filename);
    theta{i} = table2array(T(:,1));
    theta_star{i} = table2array(T(:,2));
  
    
    i=i+1;
end
    

nobus = 79
% voltage difference
figure(1)
hold on
plot(1:79, v{1}-v_star{1},'b','linestyle','-','linewidth',3);
% plot(1:79, v_star{1},'b','linestyle',':','linewidth',3);
plot(1:79, v{2}-v_star{2},'r','linestyle','-','linewidth',3);
% plot(1:79, v_star{2},'r','linestyle',':','linewidth',3);
plot(1:79, v{3}-v_star{3},'g','linestyle','-','linewidth',3);
% plot(1:79, v_star{3},'g','linestyle',':','linewidth',3);
plot(1:79, v{4}-v_star{4},'k','linestyle','-','linewidth',3);
plot(1:79, v{5}-v_star{5},'m','linestyle','-','linewidth',3);
% plot(1:79, v_star{4},'y','linestyle',':','linewidth',3);
title('voltage')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('0','0.25','0.5','0.75','1')

figure(2)
hold on
plot(1:79, theta{1}-theta_star{1},'b','linestyle','-','linewidth',3);
% plot(1:79, theta_star{1},'b','linestyle',':','linewidth',3);
plot(1:79, theta{2}-theta_star{2},'r','linestyle','-','linewidth',3);
% plot(1:79, theta_star{2},'r','linestyle',':','linewidth',3);
plot(1:79, theta{3}-theta_star{3},'g','linestyle','-','linewidth',3);
% plot(1:79, theta_star{3},'g','linestyle',':','linewidth',3);
plot(1:79, theta{4}-theta_star{4},'k','linestyle','-','linewidth',3);
plot(1:79, theta{5}-theta_star{5},'m','linestyle','-','linewidth',3);
% plot(1:79, v_star{4},'y','linestyle',':','linewidth',3);
title('theta')
xlabel('bus')
ylabel('rad')
set(gca,'Ycolor','black','FontSize',20)
legend('0','0.25','0.5','0.75','1')

figure(3)
hold on
plot(1:9, p_gen{1}-p_gen_star{1},'b','linestyle','-','linewidth',3);
% plot(1:79, p_gen_star{1},'b','linestyle',':','linewidth',3);
plot(1:9, p_gen{2}-p_gen_star{2},'r','linestyle','-','linewidth',3);
% plot(1:79, p_gen_star{2},'r','linestyle',':','linewidth',3);
plot(1:9, p_gen{3}-p_gen_star{3},'g','linestyle','-','linewidth',3);
% plot(1:79, p_gen_star{3},'g','linestyle',':','linewidth',3);
plot(1:9, p_gen{4}-p_gen_star{4},'k','linestyle','-','linewidth',3);
% plot(1:79, p_gen_star{4},'y','linestyle',':','linewidth',3);
plot(1:9, p_gen{5}-p_gen_star{5},'m','linestyle','-','linewidth',3);
title('real generation')
xlabel('generator')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('0','0.25','0.5','0.75','1')

figure(4)
hold on
plot(1:9, q_gen{1}-q_gen_star{1},'b','linestyle','-','linewidth',3);
% plot(1:79, q_gen_star{1},'b','linestyle',':','linewidth',3);
plot(1:9, q_gen{2}-q_gen_star{2},'r','linestyle','-','linewidth',3);
% plot(1:79, q_gen_star{2},'r','linestyle',':','linewidth',3);
plot(1:9, q_gen{3}-q_gen_star{3},'g','linestyle','-','linewidth',3);
% plot(1:79, q_gen_star{3},'g','linestyle',':','linewidth',3);
plot(1:9, q_gen{4}-q_gen_star{4},'k','linestyle','-','linewidth',3);
plot(1:9, q_gen{5}-q_gen_star{5},'m','linestyle','-','linewidth',3);
% plot(1:79, q_gen_star{4},'y','linestyle',':','linewidth',3);
title('reactive generation')
xlabel('generator')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('0','0.25','0.5','0.75','1')


% figure(2)
% hold on
% plot(1:79, theta,'b','linestyle','--','linewidth',3);
% plot(1:79, theta_star,'r','linestyle',':','linewidth',3);
% title('THETA')
% xlabel('bus')
% ylabel('p.u.')
% set(gca,'Ycolor','black','FontSize',20)
% legend('GML result','PowerModel')
% 
% figure(3)
% hold on
% plot(1:79, P_gen,'b','linestyle','--','linewidth',3);
% plot(1:79, P_gen_star,'r','linestyle',':','linewidth',3);
% title('Real Gen')
% xlabel('bus')
% ylabel('p.u.')
% set(gca,'Ycolor','black','FontSize',20)
% legend('GML result','PowerModel')
% 
% figure(4)
% hold on
% plot(1:79, Q_gen,'b','linestyle','--','linewidth',3);
% plot(1:79, Q_gen_star,'r','linestyle',':','linewidth',3);
% title('Reactive Gen')
% xlabel('bus')
% ylabel('p.u.')
% set(gca,'Ycolor','black','FontSize',20)
% legend('GML result','PowerModel')