clc;clear all;close all

% In this code, we apply the proposed (PGL) method for the task of joint
% data imputation and product graph learning.

% Here, we use the PM 2.5 data collected across 30 air quality monitoring
% stations in India over the year 2018. The data is available at 
% https://app.cpcbccr.com.

% We, preprocess the data and save it in the usable format in the file 
% AQI_Data_WithCoords.mat. Here, we also give the latitude and longitude
% locations as the location coordinates for ploting the figures.


addpath ./misc/
addpath ./matrix_completion/
load('AQI_Data_WithCoords.mat')

% Data preprocessing step

% We accumulate the data over T months. The data should be of dimension:
% Number of sensors times Number of months
Sq = 12; % 12 months
D = [];
for i = 0:Sq-1
    idx = i*30+1:i*30+30;
    I(i+1,:) = idx;
    Data{i+1} = X(:,idx);
    D = [D sum(Data{i+1},2)];
end
clear i I 
%%
D = D/30; % 30 days a month
X_Original = D;
X_data = D;
%% Obtain the covariance matrices from the Data
Sp = cov(X_data'); % Cov matrix for sensors
Sq = cov(X_data); % Cov matrix for time
Sp = 0.005*Sp/size(Sp,1);
Sq = 0.051*Sq/size(Sq,1);

p = size(Sp,1);
q = size(Sq,1);

param.b1 = 0.1; % ap*norm(Lp,'fro')^2
param.bp = 0.005;   % bp*trace(S*Lp)
% param.gp = 0;   %gp*trace(Vp'*Lp*Vp)

param.b2 = 1.5;   % aq*norm(Lq,'fro')^2
param.bq = 0.051;   % bq*trace(T*Lq)
% param.gq = 0;  % gq*trace(Vq'*Lq*Vq)

Sp = param.bp*Sp;
Sq = param.bq*Sq;

max_iter = 10;

% Initialization
k1 = 1; k2 = 1;
Gd = gsp_sensor(p);
Gp_i = Gd;
Gp_i = gsp_estimate_lmax(Gp_i);
[V, ~] = eig(full(Gd.L));
Vp = V(:,1:k1);
clear V Gd

Gd = gsp_sensor(q);
Gq_i = Gd;
Gq_i = gsp_estimate_lmax(Gq_i);
[V, ~] = eig(full(Gd.L));
Vq = V(:,1:k2);
clear V Gd
De_p = [];
De_q = [];
for i = 1:max_iter
    ratio = [1 0 0];
    
    prob_params.gamma_n = .1;
    prob_params.gamma_r = param.bp;
    prob_params.gamma_c = param.bq;
    solver_params.rho_ADMM = .09;
    SeedVal = 0;
    [X_data, stat_MC_KNN] = MC_function(X_data, Gq_i, Gp_i, ratio, SeedVal, prob_params, solver_params);
    X_data = double(X_data);
    
    Sp = cov(X_data'); % Cov matrix for sensors
    Sq = cov(X_data); % Cov matrix for time
    Sp = Sp/size(Sp,1);
    Sq = Sq/size(Sq,1);
    Sp = param.bp*Sp;
    Sq = param.bq*Sq;

    
   
    [Lp_i, Lq_i] = Learn_PGL(Sp, Sq, param);
    
    Wp_i = diag(diag(Lp_i)) - Lp_i;
    Wq_i = diag(diag(Lq_i)) - Lq_i;
    Wp_i(abs(Wp_i)<0.001) = 0;
    Wq_i(abs(Wq_i)<0.001) = 0;
    
    Lp_i = diag(sum(Wp_i)) - Wp_i;
    Lq_i = diag(sum(Wq_i)) - Wq_i;
    
    [V, D] = eig(Lp_i);
    Vp = V(:,1:k1);
    De_p = [De_p diag(D)];
    clear V D
    [V, D] = eig(Lq_i);
    Vq = V(:,1:k2);
    De_q = [De_q diag(D)];
    clear V D
    
    
    if i>1
        objec(i) = norm(Lp_i - Lp_prev,'fro')/norm(Lp_prev,'fro')...
            + norm(Lq_i-Lq_prev,'fro')/norm(Lq_prev);
    end
    if i>1 & objec(i)<0.0001
        break
    end
    
    Lp_prev = Lp_i;
    Lq_prev = Lq_i;
end

%%
C = coords;

coords(1,2) = 16;
coords(10,2) = 12;

Gp = gsp_graph(Wp_i);
Gp.coords = coords;
G_m = graph(Wp_i);
figure(2);clf;fp = plot(G_m);
fp.XData = coords(:,1);
fp.YData = coords(:,2);


figure(1);clf;gsp_plot_graph(Gp);
gsp_plot_signal(Gp, X_data(:,12));
names{8} = '';
names{12} = '';
names{24} = ''
names{13} = '';
names{22} = '';
names{16} = '';
names{7}= '';
names{8} = '';
names{28} = '';
names{23} = '';
names{2} = '';
names{6} = '';
names{19} = '';
names{5} = '';
names{4} = '';
names{14} = '';
text(coords(:,1), coords(:,2), names, 'VerticalAlignment','bottom','HorizontalAlignment','center','Position',[0.1,0.1,0])
Gq = graph(Wq_i);
figure(2);clf;plot(Gq,'Layout','circle')

%%
Wq_i(abs(Wq_i)<0.02) = 0;
G2 = graph(Wq_i);
figure(4);fq = plot(G2,'Layout','force');


Months = {'January','February','March','April','May','June','July','August','September','October','November','December'};
fq.NodeFontSize = 10;
fq.NodeLabel = Months;
fq.MarkerSize = 7;
G2.Nodes.value = X_data(8,:)';
G2.Nodes.NodeColors = G2.Nodes.value;
fq.NodeCData = G2.Nodes.NodeColors;
axis off
%%
I = imread('IndiaMap_2.jpg');
figure(3);imagesc(I)

rmpath ./misc/
rmpath ./matrix_completion/


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

