

% This script should give me the paper ready figures for Fig. 2 
% We will plot for the Fmeasure (Community graphs) for vaious methods
clc;clear all;close all

% Figure 1 of this code is for F-score Lp
% Figure 2 of this code is for F-score Lq
% Figure 3 of this code is for F-score Ln

% PGL Graphs
load('Result_PGL_Community_Re.mat')
figure(1);clf; semilogx(NumSignals, smooth(sum(f_p,2)/nIter), 'k','LineWidth',2); hold on
% title('Graph Lp');hold on
figure(2);clf;semilogx(NumSignals,smooth(sum(f_q,2)/nIter),'k','LineWidth',2); hold on
% title('Graph Lq');hold on
figure(3);clf;semilogx(NumSignals,smooth(sum(f,2)/nIter),'k','LineWidth',2); hold on
% title('Graph L');hold on;
figure(1);legend('PGL')
figure(2);legend('PGL')
figure(3);legend('PGL')

%%
% Dongs method
clc;clear all
load('Result_Dong_Community_Re.mat')

figure(1);semilogx(NumSignals,smooth(sum(f_p,2)/nIter),'--r', 'LineWidth', 2);hold on
figure(2);semilogx(NumSignals,smooth(sum(f_q,2)/nIter),'--r', 'LineWidth', 2); hold on
figure(3);semilogx(NumSignals,smooth(sum(f,2)/nIter),'--r', 'LineWidth', 2); hold on
figure(3);semilogx(NumSignals,smooth(sum(f_Dong,2)/nIter),'-.b','LineWidth',2); hold on

figure(1);legend('PGL','CGL')
figure(2);legend('PGL','CGL')
figure(3);legend('PGL','CGL','Dong Ln')

%%
% BigLasso method (Dotted line Purple (#bf77f6))
clc;clear all;
% load('Results_BigLasso_Community_Re.mat')
load('ResultsBiGLasso_Re_2.mat')
figure(1);semilogx(NumSignals,smooth(sum(fp,2)/nIter),':','Color','#bf77f6','LineWidth',2);hold on
figure(2);semilogx(NumSignals,smooth(sum(fq,2)/nIter),':','Color','#bf77f6','LineWidth',2); hold on
figure(3);semilogx(NumSignals,smooth(sum(f,2)/nIter),':','Color','#bf77f6','LineWidth',2); hold on
figure(1);legend('PGL','CGL','BigLASSO')
figure(2);legend('PGL','CGL','BigLASSO')
figure(3);legend('PGL','CGL','Dong Ln','BigLASSO')


%%
% BigLasso + projection (dotted lined marker X and color is purple (#bf77f6))
figure(1);semilogx(NumSignals,smooth(sum(Projfp,2)/nIter),':','Color','#bf77f6','Marker','x',...
    'MarkerSize',10,'LineWidth',2);hold on
figure(2);semilogx(NumSignals,smooth(sum(Projfq,2)/nIter),':','Color','#bf77f6','Marker','x',...
    'MarkerSize',10,'LineWidth',2);hold on
figure(3);semilogx(NumSignals,smooth(sum(Projf,2)/nIter),':','Color','#bf77f6','Marker','x',...
    'MarkerSize',10,'LineWidth',2);hold on

%%  Waheed method
clc;clear all
load('Results_Waheed_PGL.mat')

figure(1);semilogx(NumSignals,smooth(sum(f_p,2)/nIter),':', 'Color','#653700','Marker','+',...
    'MarkerSize',10, 'LineWidth', 2);hold on
figure(2);semilogx(NumSignals,smooth(sum(f_q,2)/nIter),':', 'Color','#653700','Marker','+',...
    'MarkerSize',10, 'LineWidth', 2); hold on
figure(3);semilogx(NumSignals,smooth(sum(f,2)/nIter),':', 'Color','#653700','Marker','+',...
    'MarkerSize',10, 'LineWidth', 2); hold on




%%

% figure(1);
% set( legend('$\texttt{PGL}$', '$\texttt{KronFact}$ using $\hat{\textbf{L}}_N$ from $\texttt{GL}$', '$\texttt{BiGLasso}$', '$\texttt{Projected BiGLasso}$', '$\texttt{BPGL}$', 'Location', 'Best')  ,  'Interpreter','latex');
% 
% figure(2);
% set( legend('$\texttt{PGL}$', '$\texttt{KronFact}$ using $\hat{\textbf{L}}_N$ from $\texttt{GL}$', '$\texttt{BiGLasso}$', '$\texttt{Projected BiGLasso}$','$\texttt{BPGL}$', 'Location', 'Best')  ,  'Interpreter','latex');
% 
% figure(3);
% set(legend('$\texttt{PGL}$','$\texttt{KronFact}$ using $\hat{\textbf{L}}_N$ from $\texttt{GL}$','$\textbf{L}_N$ from \texttt{GL}','$\texttt{BiGLasso}$','$\texttt{Projected BiGLasso}$','$\texttt{BPGL}$','Location', 'Best')  ,  'Interpreter','latex');
% 

figure(1);legend('PGL','CGL','BigLASSO','B+P','BPGL','Location', 'Best');set(gca,'FontSize',12);saveas(gcf,'Lp.png')
figure(2);legend('PGL','CGL','BigLASSO','B+P','BPGL','Location', 'Best');set(gca,'FontSize',12);saveas(gcf,'Lq.png')
figure(3);legend('PGL','CGL','Dong Ln','BigLASSO','B+P','BPGL','Location', 'Best');set(gca,'FontSize',12);saveas(gcf,'Ln.png')

% set( legend('$\texttt{PGL}$', '$\texttt{KronFact}$ using $\hat{\textbf{L}}_N$ from $\texttt{GL}$', '$\texttt{BiGLasso}$', '$\texttt{Projected BiGLasso}$', 'Location', 'Best')  ,  'Interpreter','latex');


% figure(1);legend('PGL','KronFact using Ln from GL','BiGLasso','Projected BigLasso', 'BPGL','Location', 'Best');set(gca,'FontSize',10);
% figure(2);legend('PGL','KronFact using Ln from GL','BiGLasso','Projected BigLasso','BPGL', 'Location', 'Best');set(gca,'FontSize',10);
% figure(3);legend('PGL','KronFact using Ln from GL','Ln from GL','BiGLasso','Projected BigLasso','BPGL','Location', 'Best');set(gca,'FontSize',10);


for i = 1:3
    
figure(i)
set(gca,'FontSize',12);
ax = gca;

xlabel('Number of samples')
ylabel('F-score')
hold on
ax.XGrid = 'on'
xlim([10 20000])
ax.Box = 'off'
% xticks([1:11])
% xticklabels({'10','50','100','250','500','1K','2K','5K', '10K', '15K', '20K'})

% legend boxoff;
end

%%
% figure(1);legend boxoff;saveas(gcf,'Lp.png');
% figure(2);legend boxoff;saveas(gcf,'Lq.png')
% figure(3);legend boxoff;saveas(gcf,'Ln.png')

% figure(1);legend boxoff;print('-depsc', 'PGL_Fm_Lp.eps');
% figure(2);legend boxoff;print('-depsc', 'PGL_Fm_Lq.eps');
% figure(3);legend boxoff;print('-depsc', 'PGL_Fm_Ln.eps');

% figure(1);print('-depsc', 'PGL_Fm_Lp.fig');
% figure(2);print('-depsc', 'PGL_Fm_Lq.fig');
% figure(3);print('-depsc', 'PGL_Fm_Ln.fig');