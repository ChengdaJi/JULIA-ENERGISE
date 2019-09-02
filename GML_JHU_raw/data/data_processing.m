clc
clear all
lf=20;
lef=16;
axf=14;
[ha, pos] = tight_subplot(2,1,[.06 .03],[.1 .01],[.1 .1]);
%% Pd data processing
Bank1_Pd=xlsread('Pd_bank1.xlsx',1,'C:G');
Bank2_Pd=xlsread('Pd_bank2.xlsx',1,'C:G');
Bank3_Pd=xlsread('Pd_bank3.xlsx',1,'C:G');
%% for time interval with 5 minutes
Pd_hour=[Bank1_Pd';Bank2_Pd';Bank3_Pd'];
size(Pd_hour)
data=[];
% increase=0:1:11;
% for h=1:24
%     data_diff{h}=(Pd_hour(:,h+1)-Pd_hour(:,h))/12;
%     data_slot{h}=Pd_hour(:,h)*ones(1,12)+data_diff{h}*increase;
%     data=[data,data_slot{h}];
% end
hour=1:25;
slots=1:1/12:25;
data=spline(hour,Pd_hour,slots);
mean_data=mean(data');
data_without=data;

data=data+0.03*diag(mean_data)*randn(15,289);


% generate Pd
Pd=[data(2:5,:);data(7:10,:);data(12:15,:)];
Pd_wo=[data_without(2:5,:);data_without(7:10,:);data_without(12:15,:)];
% generate Qd
beta=[data(1,:);data(6,:);data(11,:)];
for feeder=1:4
    Qd(feeder,:)=Pd(feeder,:).*beta(1,:);
end
for feeder=5:8
    Qd(feeder,:)=Pd(feeder,:).*beta(2,:);
end
for feeder=9:12
    Qd(feeder,:)=Pd(feeder,:).*beta(3,:);
end
Pd_1=sum(Pd(1:4,:));
Pd_2=sum(Pd(5:8,:));
Pd_3=sum(Pd(9:12,:));
Pd_1_o=sum(Pd_wo(1:4,:));
Pd_2_o=sum(Pd_wo(5:8,:));
Pd_3_o=sum(Pd_wo(9:12,:));
i=0;
for hr=0:23
    for min=0:5:55
        i=i+1;
        timeslot(i)=datenum([2015,7,1,hr,min,0]);
    end
end

subplot(ha(1))
plot(timeslot,Pd_1(1:288),'LineWidth',3,'Color','r');
hold on
plot(timeslot,Pd_2(1:288),'LineWidth',3,'Color','b');
hold on
plot(timeslot,Pd_3(1:288),'LineWidth',3,'Color','k');
hold on
plot(timeslot,Pd_1_o(1:288),'LineWidth',1,'Color','r');
hold on
plot(timeslot,Pd_2_o(1:288),'LineWidth',1,'Color','b');
hold on
plot(timeslot,Pd_3_o(1:288),'LineWidth',1,'Color','k');
grid on
datetick('x','HHPM')
ylabel('$P_{d}$ [MW]','Interpreter','latex','FontSize',lf,'Color','black');
leg=legend('Allendale bank 1','Allendale bank 2','Boolming Grove'); 
leg.FontSize=lef;

save('demand.mat','Pd','Qd');
%% Pg data processing
Pg_data1=xlsread('Pg_data.xlsx',1,'D:D')';
Pg_data2=xlsread('Pg_data.xlsx',1,'H:H')';
j=1;
for i=1:5:length(Pg_data1)
    Pg1(j)=Pg_data1(i);
    Pg2(j)=Pg_data2(i);
    j=j+1;
end
size(Pg1)
size(Pg2)
save('generation.mat','Pg1','Pg2')
subplot(ha(2))
plot(timeslot,Pg1(1:288),'-*','LineWidth',3);
hold on
plot(timeslot,Pg2(1:288),'LineWidth',3);
grid on
datetick('x','HHPM')
ylabel('$P_{solar}$ [MW/$m^{2}$]','Interpreter','latex','FontSize',lf,'Color','black');
leg=legend('Allendale','Boolming Grove'); 
leg.FontSize=lef;


%     