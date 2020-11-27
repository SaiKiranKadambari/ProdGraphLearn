function [n,S] = groupcmp(M,N)

l1 = length(M);
l2 = length(N);
M = sort(M);
N = sort(N);
n = 0;
S = zeros(0,1);

for i = 1:l1
    for j = 1:l2
        if M(i) == N(j)
            n = n+1;
            S = [S;M(i)];
        end
    end
end