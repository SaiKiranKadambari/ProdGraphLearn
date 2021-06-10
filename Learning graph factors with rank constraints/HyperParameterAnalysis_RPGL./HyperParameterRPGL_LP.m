%%%%%% Now we will compare the hyperparms analysis of RPGL, We fix the beta
%%%%%% and vary gamma. We plot the clustering accuracy varying gamma


%%%%%% This code is for LP

clc;clear all;close all
rng('default')

addpath ./misc/
load('Graphs.mat')

% Graph params 
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
N = 1000; % Fix the number of samples to 1000
nIter = 30; % number of times we want to repeat the proces
max_iter = 100; % max_iter for the proposed algorithm

param.b1 = 0.3;
param.b2 = 0.25;
param.k1 = k1;
param.k2 = k2;


gamma2 = [0.7:0.1:1.5];
gamma2 = [0.8];

for ii = 1:length(gamma1)
    for jj = 1:length(gamma2)
        
        param.g1 = gamma1(ii);
        param.g2 = gamma2(jj);
        
        for k = 1:nIter

        [X, Sp, Sq] = DataGen(Lp, Lq, N);
        Sp = Sp/N;
        Sq = Sq/N;
        [Lp_i, Lq_i, Vp, Vq, error, objec] = RPGL_Laplacian(Lp, Lq, Sp, Sq, param, max_iter);
        L_r = KronSum(Lp_i, Lq_i);
        
        
        % To evaluate the clustering performance
        [V, ~] = eigs(full(L_r) + 1e-4*eye(size(L_r)), k1*k2, 'smallestabs');
        
        [Result_Lp, est_idx_Lp] = perf_kmeans(Vp, k1, TrueIdx_Gp);
        [Result_Lq, est_idx_Lq] = perf_kmeans(Vq, k2, TrueIdx_Gq);
        [Result_L, est_idx_L] = perf_kmeans(V, k1*k2, TrueIdx_G);
        
        Lp_nmi(k) = Result_Lp(2); 
        Lq_nmi(k) = Result_Lq(2);
        L_nmi(k) = Result_L(2);
        clear Result_Lp Result_Lq est_idx_Lp est_idx_Lq
        end
        
        NMI_Lp(ii,jj) = mean(Lp_nmi);
        NMI_Lq(ii,jj) = mean(Lq_nmi);
        NMI_L(ii,jj) = mean(L_nmi);
    end
end

rmpath ./misc/
%%
zzz = smooth(NMI_Lp);
zzz = (NMI_Lp);
figure(2);clf;plot(gamma1, zzz, 'LineWidth',2); hold on
% plot(gamma1(4), zzz(4), 'rx', 'LineWidth',2, 'MarkerSize', 15)
xlabel('$\gamma_1$', 'Interpreter','latex');
ylabel('\texttt{NMI}','Interpreter','latex');
set(gca,'FontSize',11);
print('-depsc', 'HyperParam_LP.eps');


save('HyperParameterRPGL_Lp.mat')
%%

function [Lp_i, Lq_i, Vp, Vq, error, objec] = RPGL_Laplacian(Lp, Lq, Sp, Sq, param, max_iter)
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
end
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

