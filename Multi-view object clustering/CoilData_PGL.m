clc;clear all;close all

addpath ./misc/
% In this script, we need to estimate the factor graphs corresponding to
% the COIL Data without any clusters (PGL method)

load('Data_Coil.mat');
X = X';

N = size(X,2); %
p = 36; % Number of views
q = 10; % Number of objects

% Data preprocessing and Obtaining the covariance matrices
Sp = zeros(p);
Sq = zeros(q);
for i = 1:N
    x = X(:,i);
    Xi = reshape(x, p, q);
    Sp = Sp + cov(Xi');
    Sq = Sq + cov(Xi);
end
Sp = Sp/N;
Sq = Sq/N;


param.b1 = 0.51; % ap*norm(Lp,'fro')^2
param.bp = 0.15/(1000);   % bp*trace(S*Lp)
% param.gp = 0;   %gp*trace(Vp'*Lp*Vp)

param.b2 = 2;   % aq*norm(Lq,'fro')^2
param.bq = 0.0005;   % bq*trace(T*Lq)
% param.gq = 0;  % gq*trace(Vq'*Lq*Vq)

Sp = param.bp*Sp;
Sq = param.bq*Sq;

[Lp_i, Lq_i] = Learn_PGL(Sp, Sq, param);

% 
Wp = diag(diag(Lp_i)) - Lp_i;
Wp(abs(Wp)<0.1) = 0;
Gp = graph(Wp);

Wq = diag(diag(Lq_i)) - Lq_i;
Wq(abs(Wq)<0.1) = 0;
Gq = graph(Wq);

figure(1);clf;plot(Gp);title('View graph')
figure(2);clf;plot(Gq);title('Object graph')

% Save the adjacency matrices for Gephi
ZZ = Wp;
x = [1:size(ZZ,1)];
A = [x;ZZ];
y = [0 x];
A = [y' A];
% writematrix(A,'PGL_ViewGraph.csv')
clear ZZ A x y

ZZ = Wq;
x = [1:size(ZZ,1)];
A = [x;ZZ];
y = [0 x];
A = [y' A];
% writematrix(A,'PGL_ObjectGraph.csv')
clear ZZ A x y


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
[l_wf,  err] = PGL(P0, q0, C, d, 1e-6, 0.0051);
% Sanity check
% figure(100);clf;loglog(err);title('Waterfilling error')


v1 = size(Dp,2);
v2 = size(Dq,2);

Lp_i = full(reshape(Dp*l_wf(1:v1),p,p)); Lp_i(abs(Lp_i)<1e-4)=0;
Lq_i = full(reshape(Dq*l_wf(v1+1:end),q,q)); Lq_i(abs(Lq_i)<1e-4)=0;
end



function [l err] = PGL(P, q, C, d, tol, rho)

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

