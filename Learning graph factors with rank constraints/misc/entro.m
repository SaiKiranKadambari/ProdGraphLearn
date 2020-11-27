function H = entro(N,W)

H = 0;
for i = 1:length(W)
    S = length(W{i});
    H = H-S/N*log2(S/N);
end