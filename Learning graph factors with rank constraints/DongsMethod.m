clc;clear all;close all

% This code is for obtaining the graph learing performance using the method
% X.Dong,  D.Thanou,  M.Rabbat,  and  P.Frossard, “Learning graphs from
% data: A signal representation perspective,”IEEE Signal Process.Mag., 
% vol. 36, no. 3, pp. 44–63, May 2019.

% The obtained graph Laplacian matrix is product graph. To obtain the graph
% factors graph Laplacian matrices Lp and Lq from Ln is we factorize it by 
% the Rank constrained nearest Kronecker sum factorization method as in the
% paper.


addpath ./misc/

load('Graphs.mat')
L = KronSum(Lp, Lq);
W = diag(diag(L)) - L;
G = graph(W); % Product graph

[TrueIdx_Gp, binsizes] = conncomp(Gp);
[TrueIdx_Gq, binsizes] = conncomp(Gq);
[TrueIdx_G, binsizes] = conncomp(G);
clear binsizes 

p = size(Lp,1);
q = size(Lq,1);
Lp = (p/trace(Lp))*Lp;
Lq = (q/trace(Lq))*Lq;
L = KronSum(Lp, Lq);


NumSignals  = [5000]; % Maximum number of samples
% NumSignals = [10 50 100 250 500 1000 2000 5000];
nIter = 1; % number of times we want to repeat the process

% Dongs graph learning params 
param.b = 1; % sparsity param. 
param.alp = 1; % smoothness param

% KronSum Factorization params
KronParam.tp = trace(Lp); KronParam.tq = trace(Lq);
KronParam.k1 = k1; KronParam.k2 = k2;
KronParam.g1 = 30; KronParam.g2 = 30;

for i = 1:length(NumSignals)
N = NumSignals(i); % Number of graph signals
for j = 1:nIter
[X, ~, ~] = DataGen(Lp, Lq, N);
S = X*X';
S = S/N;

[L_Dong] = LearnDongGraph(S, param);

W = diag(diag(L_Dong)) - L_Dong;
W(abs(W)<0.01) = 0;
L_Dong = diag(sum(W)) - W; % L_Dong after thresholding

% Trace Normalization (May be required for some instances)
L_Dong = L_Dong/trace(L_Dong)*size(L_Dong,1); 

[V_Dong, ~] = eigs(full(L_Dong + 1e-4*eye(size(L_Dong))), k1*k2,'smallestabs');

[RLp, RLq, Vp_wf, Vq_wf, error_wf, Objec_wf] = WaterFill_RFK(L_Dong, p, q, KronParam);
Wp = diag(diag(RLp)) - RLp;
Wp(abs(Wp)<0.001) = 0;
RLp = diag(sum(Wp)) - Wp;

Wq = diag(diag(RLq)) - RLq;
Wq(abs(Wq)<0.001) = 0;
RLq = diag(sum(Wq)) - Wq;
RLp = full(RLp);
RLq = full(RLq);
L_r = KronSum(RLp, RLq);
figure(3);clf;subplot(211);loglog(error_wf);
subplot(212);loglog(Objec_wf)

% To evaluate the clustering performance
[Vp, ~] = eigs(full(RLp), k1,'smallestabs');
[Vq, ~] = eigs(full(RLq), k2,'smallestabs');
[V, ~] = eigs(full(L_r + 1e-4*eye(size(L_r))), k1*k2,'smallestabs');


[Result_Lp, est_idx_Lp] = perf_kmeans(Vp, k1, TrueIdx_Gp);
[Result_Lq, est_idx_Lq] = perf_kmeans(Vq, k2, TrueIdx_Gq);
[Result_L, est_idx_Lq] = perf_kmeans(V, k1*k2, TrueIdx_G);
[Result_L_Dong, est_idx_Lq] = perf_kmeans(V_Dong, k1*k2, TrueIdx_G);

Lp_pu(i,j) = Result_Lp(1);  Lq_pu(i,j) = Result_Lq(1); 
Lp_nmi(i,j) = Result_Lp(2); Lq_nmi(i,j) = Result_Lq(2);
Lp_ri(i,j) = Result_Lp(3);  Lq_ri(i,j) = Result_Lq(3);

L_pu(i,j) = Result_L(1);
L_nmi(i,j) = Result_L(2);
L_ri(i,j) = Result_L(3);

LDong_pu(i,j) = Result_L_Dong(1);
LDong_nmi(i,j) = Result_L_Dong(2);
LDong_ri(i,j) = Result_L_Dong(3);
clear Result_Lp Result_Lq est_idx_Lp est_idx_Lq

[~, ~, f_p(i,j), ~, ~] = graph_learning_perf_eval(Lp,RLp);
[~, ~, f_q(i,j), ~, ~] = graph_learning_perf_eval(Lq,RLq);
[~, ~, f(i,j), ~, ~] = graph_learning_perf_eval(L,L_r);
[~, ~, fDong(i,j), ~, ~] = graph_learning_perf_eval(L,L_Dong);
fprintf('Iter No (i,j) is (%d, %d) \n', i,j)
end
end
%%
Gp_rec = graph(Wp);
Gq_rec = graph(Wq);
figure(10);clf;subplot(221);plot(Gp);title('Grnd Gp')
subplot(222);plot(Gq);title('Grnd Gq')
subplot(223);plot(Gp_rec);title('Recon Gp')
subplot(224);plot(Gq_rec);title('Recon Gq')



rmpath ./misc/
%% New Functions we will use in this part of code

function [L_i] = LearnDongGraph(S, param)
% This function is the cvx based implementation of Dong formulation. 
% We will learn the product graph Ln given the covariance matrix L as input

b = param.b; % sparsity param.
alp = param.alp; % smoothness param

n = size(S,1);
v = 0.5*n*(n+1);
[A1,b1,A2,b2,M1] = laplacian_constraint_vech(n);
P_qp = b*M1'*M1;
q_qp = [alp*vec(S)'*M1]';

cvx_begin quiet
variable x(v,1)
minimize (quad_form(x,P_qp) + q_qp'*x)
subject to
A1*x == b1
A2*x <= b2
cvx_end
L_i = reshape(M1*x(1:v,1), n, n);

end