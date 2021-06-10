clc;clear all
close all

load('HyperParam.mat')


figure(1);clf;

h = surf( beta1, beta2, Result_fp');
zlabel('Fscore')
xlabel('beta1')
ylabel('beta2')
colorbar

figure(2);clf;h = surf(beta1, beta2, Result_fq');
zlabel('Fscore')
xlabel('beta1')
ylabel('beta2')
colorbar

figure(3);clf;h = surf(beta1, beta2, Result_f');
zlabel('Fscore')
xlabel('beta1')
ylabel('beta2')
colorbar


for i =1:3
    figure(i);set(gca,'FontSize',11);
end
%%
% % zlabel('F-score', 'Interpreter','latex')
% % xlabel('$\beta_1$', 'Interpreter','latex');
% % ylabel('$\beta_2$', 'Interpreter','latex');
% % xticks(beta1)
% % xticklabels(beta1)
% % yticklabels(beta2)
% colorbar
% % h.XData = beta1;
% % h.YData = beta2;
% % figure(1);set(gca,'FontSize',12);
% 
% 
% figure(2);clf;h = surf(beta2, beta1, Result_fq);
% 
% zlabel('Fscore')
% xlabel('beta1')
% ylabel('beta2')
% % zlabel('F-score', 'Interpreter','latex')
% % xlabel('$\beta_1$', 'Interpreter','latex');
% % ylabel('$\beta_2$', 'Interpreter','latex');
% % zlabel('Fscore')
% % xlabel('beta1');
% % ylabel('beta2');
% % zlabel('F score')
% % xlabel('\beta_1')
% % ylabel('\beta_2')
% % xticklabels(beta1)
% % yticklabels(beta2)
% colorbar
% % figure(2);set(gca,'FontSize',12);
% 
% 
% figure(3);clf;h = surf(beta2, beta1, Result_f);
% zlabel('Fscore')
% xlabel('beta1')
% ylabel('beta2')
% 
% % zlabel('F-score', 'Interpreter','latex')
% % xlabel('$\beta_1$', 'Interpreter','latex');
% % ylabel('$\beta_2$', 'Interpreter','latex');
% % figure(3);set(gca,'FontSize',12);
% 
% % zlabel('F-score')
% % xlabel('\beta_1')
% % ylabel('\beta_2')
% 
% % xticklabels(beta1)
% % yticklabels(beta2)
% colorbar
% 
% figure(4);clf;imagesc(beta1, beta2, Result_fp);xlabel('\beta_1');ylabel('\beta_2'); colorbar
% figure(5);clf;imagesc(beta1, beta2,Result_fq);xlabel('\beta_1');ylabel('\beta_2'); colorbar
% figure(6);clf;imagesc(beta1, beta2,Result_f);xlabel('\beta_1');ylabel('\beta_2'); colorbar
% 
% % figure(1);set(gca,'FontSize',12);
% % figure(2);set(gca,'FontSize',12);
% % figure(3);set(gca,'FontSize',12);
% % figure(1);set(gca,'FontSize',12);
% 
% for i =1:6
%     figure(i);set(gca,'FontSize',12);
% end
% 
% figure(1);print('-depsc', 'PGL_HyperParam_Surf_Lp.eps');
% figure(2);print('-depsc', 'PGL_HyperParam_Surf_Lq.eps');
% figure(3);print('-depsc', 'PGL_HyperParam_Surf_Ln.eps');
% 
% 
% figure(4);print('-depsc', 'PGL_HyperParam_ImSc_Lp.eps');
% figure(5);print('-depsc', 'PGL_HyperParam_ImSc_Lq.eps');
% figure(6);print('-depsc', 'PGL_HyperParam_ImSc_Ln.eps');