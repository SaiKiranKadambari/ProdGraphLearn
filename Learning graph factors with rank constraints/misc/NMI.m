function n = NMI(N,W,C)

L1 = length(W);
L2 = length(C);

I = 0;
for i = 1:L1
    for j = 1:L2
        [S,inter] = groupcmp(W{i},C{j});
        s1 = length(W{i});
        s2 = length(C{j});
        if s1~=0 && s2~=0 && S~=0
            I = (S/N)*log2(N*S/(s1*s2))+I;
        end
    end
end
% n = I;
n = I/((entro(N,W)+entro(N,C))/2);
% n = I/sqrt(entro(N,W)*entro(N,C));