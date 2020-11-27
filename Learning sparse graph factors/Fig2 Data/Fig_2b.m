% This script is for generating the paper ready figures for Fig. 2b
% We will plot for the Fmeasure (Community graphs) for comparison methods

clc;clear all;close all

% Figure 1 of this code is for F-score Lp
% Figure 2 of this code is for F-score Lq
% Figure 3 of this code is for F-score Ln

% PGL Graphs
load('ResultsPGL_Community.mat')
figure(1);clf; plot( smooth(sum(f_p,2)/nIter), 'k','LineWidth',2); hold on
% title('Graph Lp');hold on
figure(2);clf;plot(smooth(sum(f_q,2)/nIter),'k','LineWidth',2); hold on
% title('Graph Lq');hold on
figure(3);clf;plot(smooth(sum(f,2)/nIter),'k','LineWidth',2); hold on
% title('Graph L');hold on;
figure(1);legend('PGL')
figure(2);legend('PGL')
figure(3);legend('PGL')

%%
% Dongs method
clc;clear all
load('ResultsDong_Community.mat')

figure(1);plot(smooth(sum(f_p,2)/nIter),'r', 'LineWidth', 2);hold on
figure(2);plot(smooth(sum(f_q,2)/nIter),'r', 'LineWidth', 2); hold on
figure(3);plot(smooth(sum(f,2)/nIter),'r', 'LineWidth', 2); hold on
figure(3);plot(smooth(sum(f_Dong,2)/nIter),'-.b','LineWidth',2); hold on

figure(1);legend('PGL','CGL')
figure(2);legend('PGL','CGL')
figure(3);legend('PGL','CGL','Dong Ln')

%%
% BigLasso method (Dotted line Purple (#bf77f6))
clc;clear all;
load('ResultsBigLasso_Community.mat')
figure(1);plot(smooth(sum(fp,2)/nIter),':','Color','#bf77f6','LineWidth',2);hold on
figure(2);plot(smooth(sum(fq,2)/nIter),':','Color','#bf77f6','LineWidth',2); hold on
figure(3);plot(smooth(sum(f,2)/nIter),':','Color','#bf77f6','LineWidth',2); hold on
figure(1);legend('PGL','CGL','BigLASSO')
figure(2);legend('PGL','CGL','BigLASSO')
figure(3);legend('PGL','CGL','Dong Ln','BigLASSO')


%%
% BigLasso + projection (dotted lined marker X and color is purple (#bf77f6))
figure(1);plot(smooth(sum(Projfp,2)/nIter),':','Color','#bf77f6','Marker','x',...
    'MarkerSize',10,'LineWidth',2);hold on
figure(2);plot(smooth(sum(Projfq,2)/nIter),':','Color','#bf77f6','Marker','x',...
    'MarkerSize',10,'LineWidth',2);hold on
figure(3);plot(smooth(sum(Projf,2)/nIter),':','Color','#bf77f6','Marker','x',...
    'MarkerSize',10,'LineWidth',2);hold on

figure(1);legend('PGL','CGL','BigLASSO','B+P','Location', 'Best');set(gca,'FontSize',12);
figure(2);legend('PGL','CGL','BigLASSO','B+P','Location', 'Best');set(gca,'FontSize',12);
figure(3);legend('PGL','CGL','Dong Ln','BigLASSO','B+P','Location', 'Best');set(gca,'FontSize',12);


for i = 1:3
figure(i)
xlabel('Number of samples')
ylabel('Average F-score')
hold on
xticks([1:8])
xticklabels({'10','50','100','250','500','1000','2000','5000'})

legend boxoff;
end

figure(1);print('-depsc', 'PGL_Fm_Lp.eps');
figure(2);print('-depsc', 'PGL_Fm_Lq.eps');
figure(3);print('-depsc', 'PGL_Fm_Ln.eps');