clc;clear all;close all

load('HyperParameterRPGL_Lp.mat')

figure(1);clf;plot(gamma1, NMI_Lp, 'LineWidth', 2)
xlabel('gamma 1')
ylabel('NMI')
set(gca,'FontSize',12);


clear 
load('HyperParameterRPGL_Lq.mat')
figure(2);clf;plot(gamma2, smooth(NMI_Lq), 'LineWidth', 2)
xlabel('gamma 2')
ylabel('NMI')
set(gca,'FontSize',12);

figure(1);print('-depsc', 'HyperParam_LP.eps');
figure(2);print('-depsc', 'HyperParam_LQ.eps');