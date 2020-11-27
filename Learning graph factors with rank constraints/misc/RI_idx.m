function N = RI_idx(idx1,idx2)

L1 = max(idx1);
L2 = max(idx2);

for i = 1:L1
    W{i} = find(idx1==i);
end

for j = 1:L2
    C{j} = find(idx2==j);
end

N = RI(W,C);