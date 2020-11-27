clc;clear all; close all
% in this file we will generate the Fmeasure plots for the community graphs

% PGL Graphs
load('ResultsPGL_Community.mat')
figure(1);clf;plot(smooth(sum(f_p,2)/nIter),'LineWidth',1.5);
title('Graph Lp');hold on
figure(2);clf;plot(smooth(sum(f_q,2)/nIter),'LineWidth',1.5);
title('Graph Lq');hold on
figure(3);clf;plot(smooth(sum(f,2)/nIter),'LineWidth',1.5);
title('Graph L');hold on;
figure(1);legend('PGL')
figure(2);legend('PGL')
figure(3);legend('PGL')

% Dongs method
clc;clear all
load('ResultsDong_Community.mat')
figure(1);plot(smooth(sum(f_p,2)/nIter),'LineWidth',1.5);hold on
figure(2);plot(smooth(sum(f_q,2)/nIter),'LineWidth',1.5); hold on
figure(3);plot(smooth(sum(f,2)/nIter),'LineWidth',1.5); hold on
figure(3);plot(smooth(sum(f_Dong,2)/nIter),'LineWidth',1.5); hold on
figure(1);legend('PGL','CGL')
figure(2);legend('PGL','CGL')
figure(3);legend('PGL','CGL','Dong Ln')

% BigLasso method
clc;clear all;
load('ResultsBigLasso_Community.mat')
figure(1);plot(smooth(sum(fp,2)/nIter),'LineWidth',1.5);hold on
figure(2);plot(smooth(sum(fq,2)/nIter),'LineWidth',1.5); hold on
figure(3);plot(smooth(sum(f,2)/nIter),'LineWidth',1.5); hold on
figure(1);legend('PGL','CGL','BigLASSO')
figure(2);legend('PGL','CGL','BigLASSO')
figure(3);legend('PGL','CGL','Dong Ln','BigLASSO')

% BigLasso + Projection method
figure(1);plot(smooth(sum(Projfp,2)/nIter),'LineWidth',1.5);hold on
figure(2);plot(smooth(sum(Projfq,2)/nIter),'LineWidth',1.5); hold on
figure(3);plot(smooth(sum(Projf,2)/nIter),'LineWidth',1.5); hold on
figure(1);legend('PGL','CGL','BigLASSO','B+P')
figure(2);legend('PGL','CGL','BigLASSO','B+P')
figure(3);legend('PGL','CGL','Dong Ln','BigLASSO','B+P')




% set(gca,'FontSize',16);

