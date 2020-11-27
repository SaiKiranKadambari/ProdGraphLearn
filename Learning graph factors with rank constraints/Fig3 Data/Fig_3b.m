clc;clear all;close all

% This code is for obtaining the figure 3b of the paper.

% We give the obtained or estimated product graph factors Gp and Gq 
% form the proposed Rank-constrained product graph learning (RPGL) method.


load('ResultsFigure3.mat')

% For the Recon Lp
figure(3);clf;
ZZ = Gp_r;
hp = plot(ZZ, 'LineWidth', 5*ZZ.Edges.Weight);
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
box off
axis off

% For the Recon Lq
figure(4);clf;
ZZ = Gq_r;
hp = plot(ZZ,  'LineWidth', 7*ZZ.Edges.Weight);
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
box off
axis off

figure(3);print('-depsc', 'RPGL_Recon_Lp.eps');
figure(4);print('-depsc', 'RPGL_Recon_Lq.eps');