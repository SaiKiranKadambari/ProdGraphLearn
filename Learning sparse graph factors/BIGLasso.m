clc;clear all;close all

% In this code we will analyze the performance of the BIGLasso algorithm on
% the product graph learning task


addpath ./misc/
addpath ./BigLASSO/
load('GraphsComm.mat')
clear param Wp Wq
L = KronSum(Lp, Lq);
L = L/trace(L)*size(L,1);
W = diag(diag(L)) - L;
G = graph(W); % Product graph

figure(1);clf;subplot(211);plot(Gp);title('Gp True')
subplot(212);plot(Gq);title('Gq True')


%%

% NumSignals = [10 50 100 250 500 1000 2000 5000];
NumSignals = [50000];
nIter = 1;  % Number of iterations for one data sample

Psi0{1} = Lp; % Defining the inputs as defined in the BIGLasso code
Psi0{2} = Lq;


lambda = [0.015 0.02];
a = 60;
tol = 1e-10;
maxiter = 200;
ps(1) = size(Lp,1);
ps(2) = size(Lq,1);
tic
type = 'L1';

for i = 1:length(NumSignals)
N = NumSignals(i);
for j = 1:nIter
    [X T1 T2] = DataGen(Lp, Lq, N);
    S{1} = T1/N;
    S{2} = T2/N;
    clear T1 T2 X;

    [ Psi, Psis ] = teralasso( S, ps, type, a, tol ,lambda, maxiter);
    Z = Psi{1,1};
    Z(abs(Z)<0.001) = 0;
    Psi{1} = Z;
    clear Z
    
    Z = Psi{1,2};
    Z(abs(Z)<0.001) = 0;
    Psi{2} = Z;
    clear Z
    
    L_r = KronSum(Psi{1}, Psi{2});
    
    
    [~, ~,fp(i,j), ~, ~] = graph_learning_perf_eval(Lp, Psi{1});
    [~, ~,fq(i,j), ~, ~] = graph_learning_perf_eval(Lq, Psi{2});
    [~, ~,f(i,j), ~, ~] = graph_learning_perf_eval(L, L_r);
    
    
    % Results with the Laplacian matrix (After factorization)
    Psi_p{1} = LapProject(Psi{1});
    Psi_p{2} = LapProject(Psi{2});
    Wp = diag(diag(Psi_p{1})) - Psi_p{1};
    

    Psi_p{1}(abs(Psi_p{1})<0.1) = 0;
    Psi_p{2}(abs(Psi_p{2})<0.1) = 0;
    L_proj = KronSum(Psi_p{1}, Psi_p{2});
    
    [~, ~,Projfp(i,j), ~, ~] = graph_learning_perf_eval(Lp,Psi_p{1});
    [~, ~,Projfq(i,j), ~, ~] = graph_learning_perf_eval(Lq,Psi_p{2});
    [~, ~,Projf(i,j), ~, ~] = graph_learning_perf_eval(L, L_proj);
end
end


%%
figure(1);clf;plot(smooth(sum(fp,2)/nIter),'LineWidth',1.5);hold on
figure(2);plot(smooth(sum(fq,2)/nIter),'LineWidth',1.5); hold on
figure(3);plot(smooth(sum(f,2)/nIter),'LineWidth',1.5); hold on


% BigLasso + Projection method
figure(1);plot(smooth(sum(Projfp,2)/nIter),'LineWidth',1.5);hold on
figure(2);plot(smooth(sum(Projfq,2)/nIter),'LineWidth',1.5); hold on
figure(3);plot(smooth(sum(Projf,2)/nIter),'LineWidth',1.5); hold on
figure(1);legend('B','B+P')
figure(2);legend('B', 'B+P')
figure(3);legend('B', 'B+P')

rmpath ./misc/
rmpath ./BigLASSO/

