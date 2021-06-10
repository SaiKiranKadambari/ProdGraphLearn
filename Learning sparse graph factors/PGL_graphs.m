clc;clear all;close all
% In this coode, we will test the performance of the graph Learning methods
% in terms of Fmeasure by varing the number of available samples.

addpath ./misc/
load('GraphsComm.mat')
clear param Wp Wq
figure(1);clf;subplot(211);plot(Gp);title('Gp True')
subplot(212);plot(Gq);title('Gq True')

Lp = (P/trace(Lp))*Lp; % Trace normalization
Lq = (Q/trace(Lq))*Lq;
L = KronSum(Lp, Lq);

NumSignals  = [5000]; % Maximum number of samples
NumSignals = [10 50 100 250 500 1000 2000 5000 10000 15000 20000];
nIter = 5; % number of times we want to repeat the process

max_iter = 1; % max_iter for the proposed algorithm

param.b1 = 0.25;
param.b2 = 0.25;


for i = 1:length(NumSignals)
N = NumSignals(i); % Number of graph signals
for j = 1:nIter

[X, Sp, Sq] = DataGen(Lp, Lq, N);
Sp = Sp/N;
Sq = Sq/N;
[Lp_i, Lq_i] = Learn_PGL(Sp, Sq, param);


L_r = KronSum(Lp_i, Lq_i);

[~, ~, f_p(i,j), ~, ~] = graph_learning_perf_eval(Lp,Lp_i);
[~, ~, f_q(i,j), ~, ~] = graph_learning_perf_eval(Lq,Lq_i);
[~, ~, f(i,j), ~, ~] = graph_learning_perf_eval(L,L_r);
end
end
Wp_r = diag(diag(Lp_i)) - Lp_i;
Wq_r = diag(diag(Lq_i)) - Lq_i;
Gp_r = graph(Wp_r);
Gq_r = graph(Wq_r);


figure(2);clf;subplot(221);plot(Gp);title('Grnd Lp')
subplot(222);plot(Gq);title('Grnd Lq')
subplot(223);plot(Gp_r);title('Recon Lp');
subplot(224);plot(Gq_r);title('Recon Lq');

figure(3);clf;subplot(221);spy(Lp);title('Ground Lp');
subplot(222);spy(Lq); title('Ground Lq')
subplot(223);spy(Lp_i);title('Recon Lp')
subplot(224);spy(Lq_i);title('Recon Lq')

figure(4);clf;plot(smooth(sum(f_p,2)/nIter),'LineWidth',2); hold on
plot(smooth(sum(f_q,2)/nIter),'LineWidth', 1.5); hold on
plot(smooth(sum(f,2)/nIter),'LineWidth',1.5); hold on
title('Fmeasure')
legend('Lp','Lq','L')


rmpath ./misc/
%% Functions we will be using for this part of code
function [Lp_i, Lq_i] = Learn_PGL(Sp, Sq, param)
p = size(Sp,1);
q = size(Sq,1);

Dp = duplication_matrix(p);
Dq = duplication_matrix(q);

P0 = blkdiag(2*param.b1*Dp'*Dp, 2*param.b2*Dq'*Dq);

q1p = [vec(Sp)'*Dp]';
q1q = [vec(Sq)'*Dq]';
q0 = [q1p; q1q]; 
% constraints 
Cp = [vec(eye(p))'*Dp; kron(ones(p,1)',eye(p))*Dp];
dp = [p; zeros(p,1)];
Cq = [vec(eye(q))'*Dq; kron(ones(q,1)',eye(q))*Dq];
dq = [q; zeros(q,1)];
C = blkdiag(Cp, Cq);
d = [dp; dq];
% Solver 
[l_wf,  err] = PGL_solver(P0, q0, C, d, 1e-6, 0.0051);
% Sanity check
% figure(100);clf;loglog(err);title('Waterfilling error')


v1 = size(Dp,2);
v2 = size(Dq,2);

Lp_i = full(reshape(Dp*l_wf(1:v1),p,p)); Lp_i(abs(Lp_i)<1e-4)=0;
Lq_i = full(reshape(Dq*l_wf(v1+1:end),q,q)); Lq_i(abs(Lq_i)<1e-4)=0;
end



function [l err] = PGL_solver(P, q, C, d, tol, rho)

p = diag(P);
%mu = rand(size(C,1),1);
mu = zeros(size(C,1),1);
l_aux = p.^(-1).*(C'*mu - q);
l = max(0, l_aux);

%%
    k=1;
    res = norm(C*l - d);
    while k < 20000
        if (res > tol) %this threshold is tunable, 10^-4 is good enough
            mu = mu - rho*(C*l - d); % gradient-descent 
            l_aux = p.^(-1).*(C'*mu - q);
            l = max(0, l_aux);
        end
        k = k+1;
        res = norm(C*l - d);
        err(k) = res;
    end
end

