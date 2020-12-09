clc;clear all;close all
addpath misc\

load('AQI_Data_WithCoords.mat')
%%
% Y is the original data;
Y = X; % 360 samples of the data 
clear X

Y_true = Y; % This we will use for reference
avail_per = nnz(Y)/(size(Y,1)*size(Y,2));

p = 30; % Number of sensors
q = 12; % Number of months
T = 30; % Number of samples

for t = 1:T
    idx = t:30:360;
    Yi = Y(:, idx);
    [y_train, mask_train, y_val, mask_val, y_test, mask_test] = split_observed(Yi, [0.85, 0.15, 0], 0, 30);
    
    M_train{t,1} = mask_train;
    M_test{t,1} = mask_test;
    M_val{t,1} = mask_val;
%     r(t) = rank(D);
end
clear y_train mask_train y_val mask_val y_test mask_test idx D t 

%% Algorithm

% Use the missing data to estimate the KNN graph for initilaization

Z = [];
for i = 0:11
    idx = i*30+1:i*30+30;
    I(i+1,:) = idx;
    Data{i+1} = Y(:,idx);
    Z = [Z sum(Data{i+1},2)];
end
clear i I idx Data
%%
Z = Z/30; % 30 days a month

% Initilalization (Graphs)
Gp = gsp_nn_graph(Z);
Lp = full(Gp.L);
Lp = Lp/max(eig(Lp));

Gq = gsp_nn_graph(Z');
Lq = full(Gq.L);
Lq = Lq/max(eig(Lq));
clear Gp Gq Yi Z
%%2

param.a1 = 0.00052;
param.a2 = 0.0001;
param.b1 = 5;
param.b2 = 2;


max_iter = 20;
for k = 1:max_iter
    Sp = zeros(p,p); 
    Sq = zeros(q,q);
    Y_MC = zeros(size(Y)); 
    
    for i = 1:T
        
        idx = i:30:360; % selects the data of every month
        Yi = Y(:, idx);
        
        vec_A = double(M_train{i,1}); 
        D = diag(vec_A); % D = diag(vec(A));
        
        L = KronSum(param.a1*Lp, param.a2*Lq);
        
        B = D + L; % (D + KronSum(a1*Lp, a2*Lq))
        
        B = B + 1e-6*eye(size(B)); 
        R(i) = cond(B); 
        
        b = vec_A.*vec(Yi);
        xi = B\b; % Uses the pinv
        
        Xi = reshape(xi, p, q); % reshape to get a matrix
        
        % Reconstruction error
        y_pred = Xi(M_val{i,1});
        y_true = Yi(M_val{i,1});
        
        Nr = norm(y_pred - y_true)^2;
        Dr = norm(y_true)^2;
        Error(i) = Nr;
        Error_N(i) = Nr/Dr;
        
        RMSE(i) = sqrt(mean((y_pred - y_true).^2)); 
        % Estimate te covariance matrices for Graph Learning step
        
        Sp = Sp + cov(Xi');
        Sq = Sq + cov(Xi);
        
        Y_MC(:, idx) = Xi;
    end
    
    ranks(k,:) = R;
    T_Error(k,:) = Error;
    T_Error_N(k,:) = Error_N;
    T_RMSE(k,:) = RMSE;
    Y = Y_MC; % Update the Y for the next iteration
    
    Sp = Sp/T;
    Sq = Sq/T;
    
    Sp = param.a1*Sp;
    Sq = param.a1*Sq;
    
    [Lp, Lq] = Learn_PGL(Sp, Sq, param);
%     Lp = Lp + 10^-6*eye(p);
%     Lq = Lq + 10^-6*eye(q);
    
%     Lp = Lp/max(eig(Lp));
%     Lq = Lq/max(eig(Lq));
    
    if k>1 % Error function for convergene
        objec(k) = norm(Lp - Lp_prev,'fro')/norm(Lp_prev,'fro')...
            + norm(Lq-Lq_prev,'fro')/norm(Lq_prev)
    end
    if k>1 & objec(k)<0.0001 % Convergence
        break
    end
    
    Lp_prev = Lp;
    Lq_prev = Lq;
end


%%
I = imread('IndiaMap_2.jpg');
figure(1);clf;imagesc(I);
axis off
hold on;

Wp = diag(diag(Lp)) - Lp; 



Wp(abs(Wp)<0.05) = 0;
Gp = gsp_graph(Wp);
Gp.coords = Coords;
figure(1);
%gsp_plot_graph(Gp)
%gsp_plot_signal(Gp, X_data(:,24))
[Xg,Yg] = gplot(Gp.L, Coords);
plot(Xg,Yg,'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',8,...
    'Marker','o',...
    'Color',[0.494117647409439 0.494117647409439 0.494117647409439]);


names{8} = '';
names{12} = '';
names{24} = '';
names{13} = '';
names{22} = '';
names{16} = '';
names{7}= '';
names{8} = '';
names{28} = '';
% names{23} = '';
names{20} = '';
names{1} = '';
names{2} = '';
names{6} = '';
names{19} = '';
names{5} = '';
names{4} = '';
names{14} = '';
Txtcoords = Coords;
Txtcoords(17,2) = Coords(17,2) + 20;  % offset for overlapping text labels
Txtcoords(21,1) = Coords(21,1) + 20;  % offset for overlapping text labels
Txtcoords(3,1) = Coords(3,1) + 10;  % offset for overlapping text labels
Txtcoords(30,1) = Coords(30,1) + 10;  % offset for overlapping text labels
Txtcoords(27,1) = Coords(27,1) + 15;  % offset for overlapping text labels
text(Txtcoords(:,1), Txtcoords(:,2), names, 'VerticalAlignment','bottom',...
    'HorizontalAlignment','center','Position',[0.1,0.1,0], 'FontSize',11.5)

set(gca,'ydir','reverse')


print('-depsc', 'AQI_Data_SpaceGraph_up.eps');

%% Time graph
Wq = diag(diag(Lq)) - Lq;
Wq(abs(Wq)<0.071) = 0;
G2 = graph(Wq);
figure(2);fq = plot(G2,'Layout','force');


Months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
fq.NodeFontSize = 12;
fq.NodeLabel = Months;
fq.MarkerSize = 7;
set(gca,'Visible','off')

print('-depsc', 'AQI_Data_TimeGraph_up.eps');

z = mean(Error);
fprintf('The reconstruction error is: %0.5f \n ', z)
rmpath misc\

%% 
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
[l_wf,  err] = PGL_solver(P0, q0, C, d, 1e-6, 0.51);
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