function N = RI(W,C)

L1 = length(W);
L2 = length(C);
S = zeros(L1,L2);
inter = cell(L1,L2);

Pos = 0;
for i = 1:L1
    count(i) = length(W{i});
    Pos_tmp = count(i)*(count(i)-1)/2;
    Pos = Pos + Pos_tmp;
end

TP = 0;
for i = 1:L1
    for j = 1:L2
        [S(i,j),inter{i,j}] = groupcmp(W{i},C{j});
        if S(i,j) >= 2
            TP_tmp = S(i,j)*(S(i,j)-1)/2;
            TP = TP + TP_tmp;
        end
    end
end

FP = Pos - TP;

Neg = 0;
for i = 1:L1
    for j = 1:L1
        if i < j
            Neg_tmp = count(i)*count(j);
            Neg = Neg + Neg_tmp;
        end
    end
end

FN = zeros(1,L2);
for n = 1:L2
    for i = 1:L1
        for j = 1:L1
            if i < j
            FN_tmp(n) = S(i,n)*S(j,n);
            FN(n) = FN(n) + FN_tmp(n);
            end
        end
    end
end
FN = sum(FN);

TN = Neg - FN;

N = (TP+TN)/(TP+FP+TN+FN);