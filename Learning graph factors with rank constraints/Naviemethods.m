clc;clear all;close all
% In this code, we will test the clustering performace of the navie
% clustering approaches given the graph data. We evaluate the clustering
% performance on the kmeans and spectral clustering algotithm

% For k-means clustering algorithm,we refer to
% S. Lloyd, “Least squares quantization in pcm,” IEEE transactions on
% information theory, vol. 28, no. 2, pp. 129–137, 1982.

% For spectral clustering algorithm we refer to
% A. Y. Ng, M. I. Jordan, and Y. Weiss, “On spectral clustering: Analysisand an algorithm,” 
% in Proc. of the Adv. in neural inf. proc. sys. (NIPS),Vancouver, Canada, Dec. 2002

load('Graphs.mat')
addpath ./misc/
Lp = (p/trace(Lp))*Lp;
Lq = (q/trace(Lq))*Lq;
L = KronSum(Lp, Lq);
W = diag(diag(L)) - L;
G = graph(W); % Product graph

% [TrueIdx_Gp, binsizes] = conncomp(Gp);
% [TrueIdx_Gq, binsizes] = conncomp(Gq);
[TrueIdx_G, binsizes] = conncomp(G);
clear binsizes 

NumSignals  = [1000]; % Maximum number of samples
% NumSignals = [10 50 100 250 500 1000 2000 5000];
nIter = 1; % number of times we want to repeat the process
k = k1*k2; 

for i = 1:length(NumSignals)
N = NumSignals(i); % Number of graph signals
for j = 1:nIter
[X, ~, ~] = DataGen(Lp, Lq, N);

Xin = normalize(X, 'range', [-1 1]);

Kmeans_idx = kmeans(Xin, k, 'Replicates',30);
Kmeans_NMI(i,j) = NMI_idx(TrueIdx_G, Kmeans_idx);
Kmeans_pu(i,j) = purity_idx(TrueIdx_G, Kmeans_idx);
Kmeans_ri(i,j) = RI_idx(TrueIdx_G, Kmeans_idx);


Sc_idx = spectralcluster(Xin, k,'LaplacianNormalization', 'symmetric','ClusterMethod','kmeans');
% Sc_idx = spectralcluster(X,k);
Sc_NMI(i,j) = NMI_idx(TrueIdx_G, Sc_idx);
Sc_pu(i,j) = purity_idx(TrueIdx_G, Sc_idx);
Sc_ri(i,j) = RI_idx(TrueIdx_G, Sc_idx);

end
end

[sum(Kmeans_NMI)/nIter sum(Sc_NMI)/nIter]


rmpath ./misc/