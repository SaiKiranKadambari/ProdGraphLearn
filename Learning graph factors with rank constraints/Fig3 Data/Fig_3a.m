clc;clear all;close all

% This code is for obtaining the figure 3a of the paper.
% We give the ground truth graphs Gp and Gq that we use for 
% Product graph clustering (synthetic dataset)

load('ResultsFigure3.mat')
NodeSize = 12;

% For the Grnd thruth Lp
figure(1);clf;
ZZ = Gp;
hp = plot(ZZ, 'LineWidth', 1*ZZ.Edges.Weight);
[bins] = conncomp(ZZ);
for i = 1:4
    idx = find(bins == i);
    if i == 1;
        highlight(hp, idx, 'NodeColor', 'r');
    elseif i == 2
        highlight(hp, idx, 'NodeColor', 'b');
    elseif i == 3
        highlight(hp, idx, 'NodeColor', 'g');
    else
        highlight(hp, idx, 'NodeColor', '#7E2F8E');
    end
        
end
hp.MarkerSize = NodeSize;
hp.NodeLabel = {};
set(hp,'LineWidth',3)
box off
axis off



figure(1);print('-depsc', 'RPGL_True_Lp.eps');
figure(2);print('-depsc', 'RPGL_True_Lq.eps');
