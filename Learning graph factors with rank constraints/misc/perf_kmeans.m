function [result,est_idx] = perf_kmeans(V,K,true_idx)

%%%% evaluate the clustering performance of k-means applied to estimated 
% spectral embedding
%
% Input:
% V: representation based on which we apply k-means
% K: number of clusters
% C: groundtruth cluster indices

% Output:
% result: clustering performance
% idx: clustering index
% Measures used: purity, NMI, RI
%%%%


%% use Matlab built-in kmeans function with multiple tries = set to 20 here

[est_idx,~] = kmeans(V,K,'EmptyAction','singleton','Replicates',20);

%% compute benchmarks
pu = purity_idx(true_idx,est_idx);
nmi = NMI_idx(true_idx,est_idx);
ri = RI_idx(true_idx,est_idx);

result = [pu;nmi;ri];
