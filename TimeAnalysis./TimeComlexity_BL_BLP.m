clc;clear all;close all

% This code computes the per iteration compuatation complexity of the
% BiGLasso and Projected BiGLasso methods (with cartesian product graph 
% structure).

% A. Kalaitzis, J. Lafferty, N. D. Lawrence, and S. Zhou, “The bigraphical
% lasso,” in Proc. of the 30th Int. Conf. on Machine Learning, vol. 28,
% no. 3, Atlanta, Georgia, USA, June 2013.



addpath('misc\')

% P_var = [10:4:40];
% Q_var = 10;%[10:5:30];

Num_sig = 1000;
Nval = [10 100 1000 10000]
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
    
    
    % BigLasso
    lambda = [0.015 0.02];
    a = 60;
    tol = 1e-10;
    maxiter = 20;
    ps(1) = P;
    ps(2) = Q;
    tic
    type = 'L1';

    S{1} = Sp;
    S{2} = Sq;
    clear T1 T2;

    tic
    [ Psi, Psis ] = teralasso( S, ps, type, a, tol ,lambda, maxiter);
    t_B(k) = toc;

    tic
    Psi_p{1} = LapProject(Psi{1});
    Psi_p{2} = LapProject(Psi{2});
    t_BP(k) =  t_B(k)+toc;
    clear S ps type a tol lambda maxiter Psi Psis
end

Time_BL(1,i) = mean(t_B);
Time_BL_sc(i,:) = t_B;
Time_BLP(1,i) = mean(t_BP);
Time_BLP_sc(i,:) = t_BP;
end