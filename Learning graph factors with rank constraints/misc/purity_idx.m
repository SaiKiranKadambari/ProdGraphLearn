function p = purity_idx(idx1,idx2)

N = length(idx1);
uniqueidx1 = unique(idx1);
uniqueidx2 = unique(idx2);

% if length(uniqueidx1) ~= length(uniqueidx2)
%     error('Number of clusters does not match number of classes.');
% end;

for i = 1:length(uniqueidx1)
    W{i} = find(idx1==uniqueidx1(i));
end

for j = 1:length(uniqueidx2)
    C{j} = find(idx2==uniqueidx2(j));
end

p = purity(N,W,C);