clc;clear all;close all
% In this coode, we will test the performance of the graph Learning methods
% in terms of Fmeasure by varying the number of available samples.

% As the input graphs Gp and Gq has the commmunity structures in it, we
% should use the KronPGL (Unconstrained KronSum Factorization) to obtain
% the factor graphs


% Note.Add all the folders and subfolders of the unlocbox to the matlab
% path.


addpath ./misc/


load('GraphsComm.mat')
clear param Wp Wq
figure(1);clf;subplot(211);plot(Gp);title('Gp True')
subplot(212);plot(Gq);title('Gq True')

Lp = (P/trace(Lp))*Lp; % Trace normalization
Lq = (Q/trace(Lq))*Lq;
L = KronSum(Lp, Lq);
% L = L/trace(L)*size(L,1); 

NumSignals  = [5000]; % Maximum number of samples
% NumSignals = [10 50 100 250 500 1000 2000 5000];
nIter = 1; % number of times we want to repeat the process

param.tp = trace(Lp); % Params for KronSumFactorization
param.tq = trace(Lq);

for i = 1:length(NumSignals)
N = NumSignals(i); % Number of graph signals
for j = 1:nIter

[X, ~, ~] = DataGen(Lp, Lq, N);

Z = gsp_distanz(X').^2;
Z = Z/N;
[W_Dong] = gsp_learn_graph_l2_degrees(2.5*Z, 1.5);
W_Dong(abs(W_Dong)<0.0001) = 0;
L_Dong  = diag(sum(W_Dong)) - W_Dong;
% L_Dong = L_Dong/trace(L_Dong)*size(L_Dong,1); 

[~, ~, f_Dong(i,j), ~, ~] = graph_learning_perf_eval(L, L_Dong);

[FLp, FLq] = WaterFill_KF(L_Dong, P, Q, param);
Wp = diag(diag(FLp)) - FLp;
Wq = diag(diag(FLq)) - FLq;
Wp(abs(Wp)<0.1) = 0;
Wq(abs(Wq)<0.1) = 0;
FLp = diag(sum(Wp)) - Wp;
FLq = diag(sum(Wq)) - Wq;



FL_r = KronSum(FLp, FLq);
[~, ~, f_p(i,j), ~, ~] = graph_learning_perf_eval(Lp,FLp);
[~, ~, f_q(i,j), ~, ~] = graph_learning_perf_eval(Lq,FLq);
[~, ~, f(i,j), ~, ~] = graph_learning_perf_eval(L,FL_r);

end
end

figure(5);clf;plot(smooth(sum(f_p,2)/nIter),'LineWidth',2); hold on
plot(smooth(sum(f_q,2)/nIter),'LineWidth', 1.5); hold on
plot(smooth(sum(f,2)/nIter),'LineWidth',1.5); hold on
title('Fmeasure')
legend('Lp','Lq','L')

% Plot the Laplacian matrices
% figure(2);clf;spy(L);title('Grnd Laplacian')
% figure(3);clf;spy(L_Dong);title('Dongs Laplacian')
% figure(4);clf;plot(diag(L_Dong))


rmpath ./misc/
