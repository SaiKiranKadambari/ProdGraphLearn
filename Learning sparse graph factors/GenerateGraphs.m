clc;clear all;close all

P = 10;
Q = 15;
param.Nc = 2;
Gp = gsp_community(P, param);
param.Nc = 3;
Gq = gsp_community(Q, param);
% param.p = 1;
% param.q = 0.1;
% k1 = 2;
% k2 = 3;
% Gp = gsp_stochastic_block_graph(P, k1, param);
% Gq = gsp_stochastic_block_graph(Q, k2, param);
Lp = full(Gp.L);
Lq = full(Gq.L);

% assign edge weights to the Graphs in the range [0,1]
a = 0;
b = 1;
Wp = diag(diag(Lp)) - Lp;
ind = find(abs(Wp)>0);
r = (b-a).*rand(length(ind),1) + a;
Wp(ind) = r;
Wp = 0.5*(Wp + Wp');
Lp = diag(sum(Wp)) - Wp;
clear ind r 

Wq = diag(diag(Lq)) - Lq;
ind = find(abs(Wq)>0);
r = (b-a).*rand(length(ind),1) + a;
Wq(ind) = r;
Wq = 0.5*(Wq + Wq');
Lq = diag(sum(Wq)) - Wq;
clear a b ind r

Gp = graph(Wp);
Gq = graph(Wq);
figure(1);clf;subplot(211);plot(Gp);title('Gp')
subplot(212);plot(Gq);title('Gq')