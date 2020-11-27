clc;clear all; close all

% This code is for obtaining the graph facors Lp and Lq from the noise free
% Ln. We also report the clustering accuracy of the obtained graph factors


addpath ./misc/
load('Graphs.mat')

L = KronSum(Lp, Lq);
W = diag(diag(L)) - L;
G = graph(W); % Product graph
[TrueIdx_Gp, binsizes] = conncomp(Gp);
[TrueIdx_Gq, binsizes] = conncomp(Gq);
[TrueIdx_G, binsizes] = conncomp(G);




KronParam.tp = trace(Lp); KronParam.tq = trace(Lq);
KronParam.k1 = k1; KronParam.k2 = k2;
KronParam.g1 = 30; KronParam.g2 = 30;

[RLp, RLq, Vp_wf, Vq_wf, error_wf, Objec_wf] = WaterFill_RFK(L, p, q, KronParam);

L_r = KronSum(RLp, RLq);

% To evaluate the clustering performance
[Vp, ~] = eigs(full(RLp), k1,'smallestabs');
[Vq, ~] = eigs(full(RLq), k2,'smallestabs');
[V, ~] = eigs(full(L_r + 1e-4*eye(size(L_r))), k1*k2,'smallestabs');

[Result_Lp, est_idx_Lp] = perf_kmeans(Vp, k1, TrueIdx_Gp);
[Result_Lq, est_idx_Lq] = perf_kmeans(Vq, k2, TrueIdx_Gq);
[Result_L, est_idx_Lq] = perf_kmeans(V, k1*k2, TrueIdx_G);


i = 1; 
j = 1;
Lp_pu(i,j) = Result_Lp(1);  Lq_pu(i,j) = Result_Lq(1); 
Lp_nmi(i,j) = Result_Lp(2); Lq_nmi(i,j) = Result_Lq(2);
Lp_ri(i,j) = Result_Lp(3);  Lq_ri(i,j) = Result_Lq(3);

L_pu(i,j) = Result_L(1);
L_nmi(i,j) = Result_L(2);
L_ri(i,j) = Result_L(3);