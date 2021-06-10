clc;clear all;close all

% This code computes the per iteration compuatation complexity of the
% proposed graph factorization methods KF and RKF. 
% The results are (averaged over 20 iterations) on synthetic data as in paper.


addpath('misc\')

% P_var = [10:4:40];
% Q_var = 10;%[10:5:30];

Num_sig = 1000;
Nval = [10 100 1000 10000]
% Nval = 1000
% Nval = [100 500 1000 3000 5000 7000 10000];
% Nval = [100 200 300 400 500 1000 3000];
%%
nIter = 30;
ii = 1;
for i = 1:length(Nval)

N = Nval(i);
if N == 4
    P = 2; Q = 2;
elseif N == 10
    P = 2; Q = 5;
elseif N == 100;
    P = 10; Q = 10;
elseif N == 1000;
    P = 25; Q = 40;
else 
    P = 100; Q = 100;
end

tp = P;
tq = Q;

G = gsp_sensor(N);
L = G.L;

tic;
Ltilde = TildeTransform(L,Q,Q,P,P);
t_tilde = toc;


for k = 1:nIter
    
    tic;
    DmP = duplication_matrix(P); % NC
    DmQ = duplication_matrix(Q);

    % Constraints
    Cm1 = [vec(eye(P))'*DmP; kron(ones(P,1)',eye(P))*DmP];  % NC
    Cm2 = [vec(eye(Q))'*DmQ; kron(ones(Q,1)',eye(Q))*DmQ];  % NC
    C = [Cm1, zeros(P+1, 0.5*Q*(Q+1)); % NC
        zeros(Q+1,0.5*P*(P+1)), Cm2];

    dp = [tp ;zeros(P,1)];
    dq = [tq ;zeros(Q,1)];
    d = [tp ;zeros(P,1); tq ;zeros(Q,1)];
    Hm = 2*blkdiag( Q*DmP'*DmP,  P*DmQ'*DmQ);   
    
    
    qv = -2*[vec(eye(Q))'*Ltilde*DmP, vec(eye(P))'*Ltilde'*DmQ]';

    v =0.5*P*(P+1) + 0.5*Q*(Q+1); % Number of variable we will solve for QP
    
    [z1,  err] = PGL_solver(Hm, qv, C, d, 1e-6, 0.0051);
    
    lp_c_hat1 = z1(1:0.5*P*(P+1));
    lq_c_hat1 = z1((0.5*P*(P+1))+1:end);
    
    Lp_i = full(reshape(DmP*lp_c_hat1,P,P));  Lp_i(abs(Lp_i)<10^-4) = 0;
    Lq_i = full(reshape(DmQ*lq_c_hat1,Q,Q));  Lq_i(abs(Lq_i)<10^-4) = 0;
    t_KF(k) = toc;
    
    
    if N == 10
        tic
        [Vp, ~] = eigs(full(Lp_i)+0.001*eye(P));
        [Vq, ~] = eigs(full(Lq_i)+0.001*eye(Q));
        t_RKF(1,k) = t_KF(1,k) + toc;
    else
        tic
        [Vp, ~] = eigs(full(Lp_i)+0.001*eye(P), 5,'smallestabs');
        [Vq, ~] = eigs(full(Lq_i)+0.001*eye(Q), 5,'smallestabs');
        t_RKF(1,k) = t_KF(1,k) + toc;
    end
    
    
end


Time_KF(1,i) = mean(t_KF(10:end));
Time_KF_sc(i,:) = t_KF;
Time_RKF(1,i) = mean(t_RKF(10:end));
Time_RKF_sc(i, :) = t_RKF;

end


%%



function [l err] = PGL_solver(P, q, C, d, tol, rho)

p = diag(P);
%mu = rand(size(C,1),1);
mu = zeros(size(C,1),1);
l_aux = p.^(-1).*(C'*mu - q);
l = max(0, l_aux);

%%
k=1;
res = norm(C*l - d);
while k <= 1
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
