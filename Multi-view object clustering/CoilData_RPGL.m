clc;clear all;close all

% This code is for learning the factor graphs with k1 and k2 clusters

addpath ./misc/

load('Data_Coil.mat');
X = X';

N = size(X,2); %
p = 36; % Number of views
q = 10; % Number of objects

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
% Params for the proposed algorithm
k1 = 5; % Clusters in Objects
k2 = 3; % Clusters in Views

param.bp = 0.15/1000;   % bp*trace(S*Lp)
param.bq = 0.005;   % bq*trace(T*Lq)
Sp = param.bp*Sp;
Sq = param.bq*Sq;

param.b1 = 0.5; % ap*norm(Lp,'fro')^2
param.g1 = 5;   % gp*trace(Vp'*Lp*Vp)

param.b2 = 1;   % aq*norm(Lq,'fro')^2
param.g2 = 10;  % gq*trace(Vq'*Lq*Vq)

param.k1 = k1;
param.k2 = k2;
max_iter = 1000; % max_iter for the proposed algorithm

[Lp_i, Lq_i, Vp, Vq, error, objec] = RPGL_Laplacian(Sp, Sq, param, max_iter);
%% Ploting the results (View graph)

Wp = diag(diag(Lp_i)) - Lp_i;
Wp(abs(Wp)<0.01) = 0;
Gp = graph(Wp);
[bins,binsizes] = conncomp(Gp);
for i = 1:length(binsizes)
    id = find(bins == i)
    i
    multi_car = [];
    multi_duck = [];
    for j = 1:length(id)
        x_duck = X(id(j),:); % for Duck image
        x_car = X(72+id(j),:); % For Car Image
        I_duck = uint8(reshape(x_duck,[128, 128]));
        I_car = uint8(reshape(x_car,[128, 128]));
        multi_duck = cat(3, multi_duck, I_duck);
        multi_car = cat(3, multi_car, I_car);
    end
    h_car = figure();montage(multi_car);
    saveas(h_car,sprintf('Car_View_fig%d.png',i)); % will create FIG1, FIG2,...
    
    h_duck = figure();montage(multi_duck);
    saveas(h_duck,sprintf('Duck_View_fig%d.png',i)); % will create FIG1, FIG2,...
%     title(['Views corresponding to the cluster index : ' num2str(i)])
end

%% Ploting the results (Object graph)
Wq = diag(diag(Lq_i)) - Lq_i;
Wq(abs(Wq)<0.01) = 0;

Gq = graph(Wq);
[bins,binsizes] = conncomp(Gq);

len = 1:36:size(X,1);
for i = 1:length(binsizes)
    id = find(bins == i)
    i
    multi = [];
    for j = 1:length(id)
%         k = len
        x = X(len(id(j)),:);
        I = uint8(reshape(x,[128, 128]));
        multi = cat(3, multi, I);
    end
    h = figure();montage(multi);
    
    saveas(h,sprintf('Object_fig%d.png',i)); % will create FIG1, FIG2,...
%     title(['Objects corresponding to the cluster index : ' num2str(i)])
end

%% Save the Adjacency matrices as csv files.

% View graph
ZZ = Wp;
x = [1:size(ZZ,1)];
A = [x;ZZ];
y = [0 x];
A = [y' A];
writematrix(A,'RPGL_ViewGraph.csv')
clear ZZ A x y

% Object graph
ZZ = Wq;
x = [1:size(ZZ,1)];
A = [x;ZZ];
y = [0 x];
A = [y' A];
writematrix(A,'RPGL_ObjectGraph.csv')
clear ZZ A x y




L_r = KronSum(Lp_i, Lq_i);
ZZ = diag(diag(L_r)) - L_r;
x = [1:size(ZZ,1)];
A = [x;ZZ];
y = [0 x];
A = [y' A];
writematrix(A,'RPGL_ProductGraph.csv')
clear ZZ A x y

rmpath ./misc/




%% Functions we will use in this code
function [Lp_i, Lq_i, Vp, Vq, error, objec] = RPGL_Laplacian(Sp, Sq, param, max_iter)
p = size(Sp,1);
q = size(Sq,1);

[Lp_i, Lq_i] = Learn_PGL(Sp, Sq, param); % initial step
[Vp, ~] = eigs(full(Lp_i), param.k1,'smallestabs');
[Vq, ~] = eigs(full(Lq_i), param.k2,'smallestabs');
% Gp_i = gsp_sensor(p);
% Gq_i = gsp_sensor(q);
% [Vp, ~] = eigs(full(Gp_i.L), param.k1,'smallestabs');
% [Vq, ~] = eigs(full(Gq_i.L), param.k2,'smallestabs');
% clear Gp_i Gq_i

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
    [l_wf,  err] = PGL(P0, q0, C, d, 1e-6, 0.0051);
    v1 = size(Dp,2);
    v2 = size(Dq,2);
    
    
    Lp_i = full(reshape(Dp*l_wf(1:v1),p,p)); Lp_i(abs(Lp_i)<1e-4)=0;
    Lq_i = full(reshape(Dq*l_wf(v1+1:end),q,q)); Lq_i(abs(Lq_i)<1e-4)=0;
    
    [Vp, ~] = eigs(full(Lp_i)+1e-4*eye(size(Lp_i)), param.k1,'smallestabs');
    [Vq, ~] = eigs(full(Lq_i)+1e-4*eye(size(Lq_i)), param.k2,'smallestabs');
    
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
    
%     fprintf('trace(Vp^T*Lp_true*Vp) is %0.3f and trace(Vq^T*Lq_true*Vq) is %0.3f \n', trace(Vp'*Lp*Vp), trace(Vq'*Lq*Vq))
%     fprintf('rank(Lp) is %d and rank of (Lq) is %d \n', rank(full(Lp_i),10^-3), rank(full(Lq_i),10^-3))
    
    Lp_prev = Lp_i;
    Lq_prev = Lq_i;
end
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

