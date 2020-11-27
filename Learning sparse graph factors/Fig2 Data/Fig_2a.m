clc;clear all;close all

% This code is for generating the Fig 2a in paper

% This ResultsPGL_Community.mat file has the ground truth graphs Gp, Gq and
% reconstructed graphs.

load('ResultsPGL_Community.mat') 
% Grnd truth Gp
figure(1);clf;
ZZ = Gp;
hp = plot(ZZ, 'LineWidth', 4*ZZ.Edges.Weight);
% [bins] = conncomp(ZZ);
idx_1 = [1:4];
idx_2 = [5:10];
hp.MarkerSize = 12;
highlight(hp,idx_1,'NodeColor','b')
highlight(hp,idx_2,'NodeColor','r')
hp.NodeLabel = {};
box off
axis off


print('-depsc', 'PGLGrnd_Lp.eps');
%%
% Grnd truth Gq
figure(2);clf;
ZZ = Gq;
hp = plot(ZZ, 'LineWidth', 4*ZZ.Edges.Weight);
% [bins] = conncomp(ZZ);
idx_1 = [1:5];
idx_2 = [6:11];
idx_3 = [12:15];
hp.MarkerSize = 12;
highlight(hp,idx_1,'NodeColor','b')
highlight(hp,idx_2,'NodeColor','r')
highlight(hp,idx_3,'NodeColor','g')
hp.NodeLabel = {};
box off
axis off

print('-depsc', 'PGLGrnd_Lq.eps');

%%
% Recon Gp
figure(3);clf;
ZZ = Gp_r;
hp = plot(ZZ, 'LineWidth', 10*ZZ.Edges.Weight);
% [bins] = conncomp(ZZ);
idx_1 = [1:4];
idx_2 = [5:10];
hp.MarkerSize = 12;
highlight(hp,idx_1,'NodeColor','b')
highlight(hp,idx_2,'NodeColor','r')
hp.NodeLabel = {};
box off
axis off


print('-depsc', 'PGLRecon_Lp.eps');

%%
% Recon Gq
figure(4);clf;
ZZ = Gq_r;
hp = plot(ZZ, 'LineWidth', 9*ZZ.Edges.Weight);
% [bins] = conncomp(ZZ);
idx_1 = [1:5];
idx_2 = [6:11];
idx_3 = [12:15];
hp.MarkerSize = 12;
highlight(hp,idx_1,'NodeColor','b')
highlight(hp,idx_2,'NodeColor','r')
highlight(hp,idx_3,'NodeColor','g')
hp.NodeLabel = {};
box off
axis off

print('-depsc', 'PGLRecon_Lq.eps');


%%
figure(5);clf;spy(Lp);
figure(6);clf;spy(Lp_i)