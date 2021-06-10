clc;clear all;close all

% This code computes the per iteration compuatation complexity of the
% state-of-the-art graph learning (ignoring the product graph structre) 
% methods with and without spectral constraints.

% X. Dong, D. Thanou, M. Rabbat, and P. Frossard, “Learning graphs
% from data: A signal representation perspective,” IEEE Signal Process.
% Mag., vol. 36, no. 3, pp. 44–63, May 2019.

% S. Kumar, J. Ying, J. V. de Miranda Cardoso, and D. Palomar, “Structured
% graph learning via Laplacian spectral constraints,” in Proc. of Advs.
% in Neural Inf. Proc. Sys. (NeurIPS), Vancouver, Canada, Dec. 2019.


% The results are averaged over 20 iterations on synthetic data as in paper


addpath('misc\')

% P_var = [10:4:40];
% Q_var = 10;%[10:5:30];
Num_sig = 1000;
Nval = [10 100 1000 10000]
% Nval = [5000]
% Nval = [100 500 1000 3000 5000 7000 10000];
% Nval = [100 200 300 400 500 1000 3000];
%%
nIter = 20;
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

   
    
    param_PGL.b1 = 0.25; % param for PGL
    param_PGL.b2 = 0.25; % param for PGL
    
    for k = 1:nIter
    % Data on the graph
    X = rand(N, Num_sig);
    X = normc(X);
    Sp = 0;
    Sq = 0;
    for ii = 1:Num_sig
        x = X(:,ii);
        Xi = reshape(x, P, Q);
        Sp = Sp + cov(Xi');
        Sq = Sq + cov(Xi);
    end
    clear x Xi ii
    
    
    % Dong
    tic
    [L_wf] = GraphLearn(X, 0.1, 0.01, 1e-6, 0.0051);
    t_Dong(1,k) = toc;
    
    tic;
    if N == 4
        [Vp, ~] = eigs(full(L_wf)+0.001*eye(N), 2,'smallestabs');
    else
        [Vp, ~] = eigs(full(L_wf)+0.001*eye(N), 10,'smallestabs');
    end
    t_Pol(1,k) = t_Dong(1,k) + toc;
    end
    
%     Time_PGL(1,i) = mean(t_PGL);
%     Time_RPGL(1,i) = mean(t_RPGL);
%     Time_BL(1,i) = mean(t_B);
    Time_Dong(1,i) = mean(t_Dong);
    Time_Dong_sc(i,:) = t_Dong;
    Time_Pol(1,i) = mean(t_Pol);
    Time_Pol_sc(i,:) = t_Pol;
end


%%
figure(100);clf;
loglog(Nval, smooth(Time_Dong), 'LineWidth',2); hold on
loglog(Nval, smooth(Time_Pol), 'LineWidth',2); hold on
legend('Dong GL','Plo GL', 'Location', 'best')

ax = gca;
xlabel('Number of samples')
ylabel('F-score')
hold on
ax.YGrid = 'on'
% xlim([500 10000])
ax.Box = 'off'
legend('Dong GL','Plo GL', 'Location', 'best')

%% Functions we will be using for this part of code
function [L_wf] = GraphLearn(X, g1, g2, tol, rho)
N = size(X,1);
Dm = sparse(duplication_matrix(N));
% lp_c = abs(vech(Lp)); % vec(L) = Dm*abs(lp_c);
% figure(3); spy(reshape(Dm*lp_c,N,N))
S = g1*X*X';
beta = g2;

P = 2*beta*(Dm'*Dm);
q =(vec(S)'*Dm)';
C = [vec(eye(N))'*Dm; kron(sparse(ones(N,1))',speye(N))*Dm];
d = [N; zeros(N,1)];

[l_wf err] = PGL_solver(P, q, C, d, tol,rho);

L_wf = reshape(Dm*l_wf,N,N); L_wf(abs(L_wf)<10^(-3))=0;
% figure(6); spy(L_wf);
% figure(7); loglog(err);
end



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
