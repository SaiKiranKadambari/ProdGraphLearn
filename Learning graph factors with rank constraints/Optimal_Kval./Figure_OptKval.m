clc;clear all;close all


load('OptimalK_val_1.mat')

zzz_Q = mean(Objec_LQ,2);
figure(1);clf;plot(kVal,zzz_Q,'k', 'LineWidth', 2); hold on
plot(kVal(3),zzz_Q(3),'rx', 'MarkerSize', 15, 'LineWidth', 2);hold on
xlabel('Kp')
ylabel('Objective fun.')
set(gca,'FontSize',13);
figure(1);print('-depsc', 'Optimal_KP_Objec.eps');

zzz2 = mean(Res_Lq,2);
figure(2);clf;plot(kVal,(zzz2),'k', 'LineWidth', 2); hold on
plot(kVal(3),zzz2(3),'rx', 'MarkerSize', 15, 'LineWidth', 2); hold on 
xlabel('Kp')
ylabel('Trace Lp')
set(gca,'FontSize',13);
figure(2);print('-depsc', 'Optimal_KP_Res.eps');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
load('OptimalK_val_2.mat')

clc
zzz_Q = (mean(Objec_LQ,2));
figure(3);plot(kVal,(zzz_Q),'k', 'LineWidth', 2); hold on
plot(kVal(4),zzz_Q(4),'rx', 'MarkerSize', 15, 'LineWidth', 2);hold on
xlabel('KQ')
ylabel('Objective fun.')
set(gca,'FontSize',13);
figure(3);print('-depsc', 'Optimal_KQ_Objec.eps');


zzz2 = mean(Res_Lp,2);
figure(4);plot(kVal,(zzz2),'k', 'LineWidth', 2); hold on
plot(kVal(4),zzz2(4),'rx', 'MarkerSize', 15, 'LineWidth', 2); hold on 
xlabel('KQ')
ylabel('Trace')
set(gca,'FontSize',13);
figure(4);print('-depsc', 'Optimal_KQ_Res.eps');