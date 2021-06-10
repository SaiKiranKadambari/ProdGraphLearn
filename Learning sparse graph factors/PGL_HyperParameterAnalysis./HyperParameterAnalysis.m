clc;clear all;close all


% In this code the aim is to analyze the hyperparmeter analysis of the PGL
% method on synthetic data. We fix the number of available training samples
% and vary beta1 and beta 2. We record f-measure averaged over 15 
% iterations and plot them



addpath ./misc/
load('GraphsComm.mat')
clear param Wp Wq
% figure(1);clf;subplot(211);plot(Gp);title('Gp True')
% subplot(212);plot(Gq);title('Gq True')

Lp = (P/trace(Lp))*Lp; % Trace normalization
Lq = (Q/trace(Lq))*Lq;
L = KronSum(Lp, Lq);

N = 500; % Fix the number of samples to 500
nIter = 15; % number of times we want to repeat the proces

beta1 = [0.1:0.02:0.5];
beta2 = [0.1:0.02:0.4];
% beta1 = [0.1 0.25 0.05];
% beta2 = [0.1 0.25];



    
for ii = 1:length(beta1)
param.b1 = beta1(ii);
for jj = 1:length(beta2)
    [ii jj]
param.b2 = beta2(jj);


for j = 1:nIter

    [X, Sp, Sq] = DataGen(Lp, Lq, N);
    Sp = Sp/N;
    Sq = Sq/N;
    [Lp_i, Lq_i] = Learn_PGL(Sp, Sq, param);


    L_r = KronSum(Lp_i, Lq_i);

    [~, ~, f_p(j), ~, ~] = graph_learning_perf_eval(Lp,Lp_i);
    [~, ~, f_q(j), ~, ~] = graph_learning_perf_eval(Lq,Lq_i);
    [~, ~, f(j), ~, ~] = graph_learning_perf_eval(L,L_r);
end

Result_fp(ii,jj) = mean(f_p);
Result_fq(ii,jj) = mean(f_q);
Result_f(ii,jj) = mean(f);


% Result_p(:,:,j) = {f_p};
% Result_q(:,:,j) = {f_q};
% Result_n(:,:,j) = {f};


end
end


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

