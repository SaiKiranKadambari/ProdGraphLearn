function p = purity(N,W,C)

L1 = length(W);
L2 = length(C);

for i = 1:L1
    tmp = 0;
    for j = 1:L2
        tmp = max(groupcmp(W{i},C{j}),tmp);
    end
    S(i) = tmp;
end

p = sum(S)/N;