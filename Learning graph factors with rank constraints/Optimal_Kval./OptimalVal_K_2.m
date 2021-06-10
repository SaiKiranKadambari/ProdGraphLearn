%%%% Here we vary K2 and try to find the ootimizal value. The best value
%%%% sjhould be 4


clc;clear all;close all
rng('default')
addpath ./misc/
load('Graphs.mat')
p = size(Lp,1);
q = size(Lq,1);
Lp = (p/trace(Lp))*Lp;
Lq = (q/trace(Lq))*Lq;
L = KronSum(Lp, Lq);
W = diag(diag(L)) - L;
G = graph(W); % Product graph

[TrueIdx_Gp, binsizes] = conncomp(Gp);
[TrueIdx_Gq, binsizes] = conncomp(Gq);
[TrueIdx_G, binsizes] = conncomp(G);
clear binsizes W Wp Wq
%%

kVal = [1:7];
NumSignals  = [1000]; % Maximum number of samples
% NumSignals = [10 50 100 250 500 1000 2000 5000];
nIter = 50; % number of times we want to repeat the process
param.max_iter = 20000; % max number of iterations for waterfilling solver.
max_iter = 100; % max_iter for the proposed algorithm
param.b1 = 0.25;
param.g1 = 0.2;
param.k1 = k1;

param.b2 = 0.25;
param.g2 = 0.5;
param.k2 = k2;


for ii = 1:length(kVal)
    for j = 1:nIter
        param.k1 = kVal(ii);
%         param.k2 = kVal(ii);
        [X, Sp, Sq] = DataGen(Lp, Lq, NumSignals);
        Sp = Sp/NumSignals;
        Sq = Sq/NumSignals;

        [Lp_i, Lq_i, Vp, Vq, error, objec, Residual_Lp, Residual_Lq] = RPGL_Laplacian(Lp, Lq, Sp, Sq, param, max_iter);
        Objec_LQ(ii,j) = objec(length(objec));
        Res_Lp(ii,j) = Residual_Lp(length(Residual_Lp));
%         Res_Lq(ii,j) = Residual_Lq(length(Residual_Lq));
%         clear objec Residual_Lp Residual_Lq
    end   
end

rmpath ./misc/
%%
clc
zzz_Q = (mean(Objec_LQ,2));
figure(1);clf;plot(kVal,(zzz_Q),'k', 'LineWidth', 2); hold on
plot(kVal(4),zzz_Q(4),'rx', 'MarkerSize', 15, 'LineWidth', 2);hold on

% title('Objective fun. varying $K_P$',  'Interpreter','latex');
xlabel('$K_Q$', 'Interpreter','latex')
ylabel('\textbf{Objective fun. \texttt{RPGL}}',  'Interpreter','latex')
set(gca,'FontSize',13);
%%
zzz2 = mean(Res_Lp,2);
figure(2);clf;plot(kVal,(zzz2),'k', 'LineWidth', 2); hold on
plot(kVal(4),zzz2(4),'rx', 'MarkerSize', 15, 'LineWidth', 2); hold on 

% title('${\mathrm{tr}}(\textbf{V}^T_P \,\,  \textbf{L}_P \textbf{V}_P )$',  'Interpreter','latex');
xlabel('$K_Q$', 'Interpreter','latex')
ylabel('${\mathrm{tr}}(\textbf{V}^T_Q \textbf{L}_Q \textbf{V}_Q )$',  'Interpreter','latex')
set(gca,'FontSize',13);

%%
figure(1);print('-depsc', 'Optimal_KQ_Objec.eps');
figure(2);print('-depsc', 'Optimal_KQ_Res.eps');

%%
function [Lp_i, Lq_i, Vp, Vq, error, objec, Residual_Lp, Residual_Lq] = RPGL_Laplacian(Lp, Lq, Sp, Sq, param, max_iter)
p = size(Sp,1);
q = size(Sq,1);

Gp_i = gsp_sensor(p);
Gq_i = gsp_sensor(q);
[Vp, ~] = eigs(full(Gp_i.L), param.k1,'smallestabs');
[Vq, ~] = eigs(full(Gq_i.L), param.k2,'smallestabs');
clear Gp_i Gq_i

Dp = duplication_matrix(p);
Dq = duplication_matrix(q);

P0 = blkdiag(2*param.b1*Dp'*Dp, 2*param.b2*Dq'*Dq);
q1p = [vec(Sp)'*Dp]';
q1q = [vec(Sq)'*Dq]';
% q2p = [param.g1*vec(Vp*Vp')'*Dp]';
% q2q = [param.g2*vec(Vq*Vq')'*Dq]';
% q0 = [q1p + q2p; q1q + q2q];

% Constraints
Cp = [vec(eye(p))'*Dp; kron(ones(p,1)',eye(p))*Dp];
dp = [p; zeros(p,1)];
Cq = [vec(eye(q))'*Dq; kron(ones(q,1)',eye(q))*Dq];
dq = [q; zeros(q,1)];
C = blkdiag(Cp, Cq);
d = [dp; dq];

clear Cp Cq dp dq

q1p = [vec(Sp)'*Dp]';
q1q = [vec(Sq)'*Dq]';
Ep = [];
Eq = [];
for k = 1:max_iter
    
    q2p = [param.g1*vec(Vp*Vp')'*Dp]';
    q2q = [param.g2*vec(Vq*Vq')'*Dq]';
    q0 = [q1p + q2p; q1q + q2q];
    
%     [l_wf err] = waterfill_solver(P0, q0, C, d, param.max_iter);
    [l_wf,  err] = PGL_solver(P0, q0, C, d, 1e-6, 0.0051);
    v1 = size(Dp,2);
    v2 = size(Dq,2);
    
    
    Lp_i = full(reshape(Dp*l_wf(1:v1),p,p)); Lp_i(abs(Lp_i)<1e-4)=0;
    Lq_i = full(reshape(Dq*l_wf(v1+1:end),q,q)); Lq_i(abs(Lq_i)<1e-4)=0;
    
    [Vp, ~] = eigs(full(Lp_i), param.k1,'smallestabs');
    [Vq, ~] = eigs(full(Lq_i), param.k2,'smallestabs');
    
    Ep = [Ep eig(full(Lp_i))];
    Eq = [Eq eig(full(Lq_i))];
    
    % convergence
    if k>1
        error(k) = norm(Lp_prev- Lp_i,'fro') + norm(Lq_prev-Lq_i,'fro');
        objec(k) = trace(Sp*Lp_i) + trace(Sq*Lq_i) ...
            + param.b1*norm(Lp_i,'fro')^2 + param.b2*norm(Lq_i,'fro')^2 ...
            + param.g1*trace(Vp'*Lp_i*Vp) + param.g2*trace(Vq'*Lq_i*Vq);
    end
    if k>1 & error(k)<1e-4
        fprintf('Converged at %d iteration \n',k)
        break;
    end
    
    fprintf('trace(Vp^T*Lp_true*Vp) is %0.3f and trace(Vq^T*Lq_true*Vq) is %0.3f \n', trace(Vp'*Lp*Vp), trace(Vq'*Lq*Vq))
    fprintf('rank(Lp) is %d and rank of (Lq) is %d \n', rank(full(Lp_i),10^-3), rank(full(Lq_i),10^-3))
    
    Lp_prev = Lp_i;
    Lq_prev = Lq_i;
    
    Residual_Lp(k) = trace(Vp'*Lp_i*Vp);
    Residual_Lq(k) = trace(Vq'*Lq_i*Vq);
end
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

