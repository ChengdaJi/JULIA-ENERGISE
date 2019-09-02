clc
clear all
% alpha
alpha=xlsread('alpha.xlsx');
Lambda=xlsread('Lambda_West_20180731.xlsx');
alpha_10=xlsread('alpha_10.xlsx');
save('price.mat','alpha','Lambda','alpha_10');