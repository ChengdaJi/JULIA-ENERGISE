% Network Data
Pd_net=[0.336608, 1.14398, 1.32062, 0.535618];
Qd_net = [0.208629, 0.708931, 0.818389, 0.331949];

% PowerModels result
P_gen_star = [1.23518, 0, 0, 2.12222];
Q_gen_star = [0.968725, 0, 0, 0.82532];
voltage_star = [1, 0.981131, 0.975857, 1];
theta_star = [0, -0.00993128, -0.0204045, 0.0213993];

% GML result
P_gen = [1.225447093, 0, 0, 2.122222204];
Q_gen = [2.396199979, 0, 0, 1.794304374];
voltage = [1, 0.928504872, 1.040862763, 1];
theta = [0, -0.021402648, -0.005909534, 0.02139931];

% voltage difference
figure(1)
hold on
plot(1:4, voltage,'b','linestyle','-','linewidth',3);
plot(1:4, voltage_star,'r','linestyle',':','linewidth',3);
title('voltage')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('GML result','PowerModel')


figure(2)
hold on
plot(1:4, theta,'b','linestyle','--','linewidth',3);
plot(1:4, theta_star,'r','linestyle',':','linewidth',3);
title('THETA')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('GML result','PowerModel')

figure(3)
hold on
plot(1:4, P_gen,'b','linestyle','--','linewidth',3);
plot(1:4, P_gen_star,'r','linestyle',':','linewidth',3);
title('Real Gen')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('GML result','PowerModel')

figure(4)
hold on
plot(1:4, Q_gen,'b','linestyle','--','linewidth',3);
plot(1:4, Q_gen_star,'r','linestyle',':','linewidth',3);
title('Reactive Gen')
xlabel('bus')
ylabel('p.u.')
set(gca,'Ycolor','black','FontSize',20)
legend('GML result','PowerModel')


% P_inj =  -Pd_net + P_gen;
% Q_inj = -Qd_net + Q_gen;
% 
% [diff,diffreactive] = nodal_balance_lin(P_inj, Q_inj, theta, voltage,...
% G,B)
% 
% 
% 
% function [diff_real,diff_reactive] = nodal_balance(P,Q,theta,v,G,B)
% 
% diff_real = zeros(1,4);
% diff_reactive = zeros(1,4);
% loss_real = zeros(1,4);
% loss_reactive = zeros(1,4);
% for bus = 1:4
%     for bus_other = 1:4
%         loss_real(bus) = v(bus)*v(bus_other)*...
%             (G(bus,bus_other)*cos(theta(bus)-theta(bus_other))...
%             +B(bus,bus_other)*sin(theta(bus)-theta(bus_other)));
%         loss_reactive(bus) = v(bus)*v(bus_other)*...
%             (G(bus,bus_other)*sin(theta(bus)-theta(bus_other))...
%             -B(bus,bus_other)*cos(theta(bus)-theta(bus_other)));
%     end
%     diff_real(bus) = -P(bus)+loss_real(bus);
%     diff_reactive(bus) = -Q(bus)+loss_reactive(bus);
% end
% 
% 
% end
% 
% function [diff_real,diff_reactive] = nodal_balance_lin(P,Q,theta,v,G,B)
% 
% Pd_net=[0.336608, 1.14398, 1.32062, 0.535618];
% Qd_net = [0.208629, 0.708931, 0.818389, 0.331949];
% P_gen_star = [1.23518, 0, 0, 2.12222];
% Q_gen_star = [0.968725, 0, 0, 0.82532];
% v_star = [1, 0.981131, 0.975857, 1];
% theta_star = [0, -0.00993128, -0.0204045, 0.0213993];
% P_star =  -Pd_net + P_gen_star;
% Q_star = -Qd_net + Q_gen_star;
% 
% diff_real = zeros(1,4);
% diff_reactive = zeros(1,4);
% loss_real = zeros(1,4);
% loss_reactive = zeros(1,4);
% for bus = 1:4
%     for bus_other = 1:4
%         loss_real(bus) = (v(bus)-v_star(bus))*v_star(bus_other)*...
%             (G(bus,bus_other)*cos(theta_star(bus)-theta_star(bus_other))...
%             +B(bus,bus_other)*sin(theta_star(bus)-theta_star(bus_other)))...
%             +(v(bus_other)-v_star(bus_other))*v(bus)*...
%             (G(bus,bus_other)*cos(theta_star(bus)-theta_star(bus_other))...
%             +B(bus,bus_other)*sin(theta_star(bus)-theta_star(bus_other)))...
%             +v_star(bus)*v_star(bus_other)*(theta(bus)-theta_star(bus))*...
%             (-G(bus,bus_other)*cos(theta_star(bus)-theta_star(bus_other))...
%             +B(bus,bus_other)*sin(theta_star(bus)-theta_star(bus_other)))...
%             +v_star(bus)*v_star(bus_other)*(theta(bus_other)-theta_star(bus_other))*...
%             (G(bus,bus_other)*cos(theta_star(bus)-theta_star(bus_other))...
%             -B(bus,bus_other)*sin(theta_star(bus)-theta_star(bus_other)));
%         loss_reactive(bus) = v(bus)*v(bus_other)*...
%             (G(bus,bus_other)*sin(theta(bus)-theta(bus_other))...
%             -B(bus,bus_other)*cos(theta(bus)-theta(bus_other)));
%     end
%     diff_real(bus) = -P(bus)+P_star(bus)+loss_real(bus);
%     diff_reactive(bus) = -Q(bus)+loss_reactive(bus);
% end
% 
% 
% end
