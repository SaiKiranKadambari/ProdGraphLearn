clc;clear all;close all

% This code is for generating the figure 3c of the paper.
% In Fig. 3c, we show the convergence of error of RPGL and R-KronFact which use alternating minimization. 





%% RPGL - given data(learning part) 1000 iters plot error (dashed black)
load('ResultsRPGL_N1000.mat')
figure(3);clf;loglog(error(2:6),'--k','LineWidth',2);hold on
xlabel('Number of iterations')
ylabel('Error')
set(gca,'FontSize',12);
% print('-depsc', 'RPGL_ErrorPlot.eps');
% RPGL - given data(learning part) 1000 iters plot error (dashed black)
% figure(4);clf;loglog(objec(2:end),'--k','LineWidth',2);
% xlabel('Number of iterations')
% ylabel('Objective function')
% set(gca,'FontSize',12);
% print('-depsc', 'RPGL_ObjecPlot.eps');

%%  RPGL method run 1000 iterations and plot the error (Red dashed marker + ).
clear all;clc
load('Results_KronFact.mat')
figure(3);loglog(error_wf(2:160), '--','Color','r','Marker','+',...
    'MarkerSize',10,'LineWidth',2);


xlabel('Number of iterations')
ylabel('Error')
legend('RPGL','RPFL with Dong')
legend boxoff;

set(gca,'FontSize',12);
print('-depsc', 'ErrorPlot.eps');
% 
% figure(2);clf;loglog(Objec_wf(2:end), '--','Color','r','Marker','+',...
%     'MarkerSize',10,'LineWidth',2)
% xlabel('Number of iterations')
% ylabel('Objective function')
% set(gca,'FontSize',12);
% print('-depsc', 'KronFactor_ObjecPlot.eps');