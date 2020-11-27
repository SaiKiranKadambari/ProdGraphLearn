function [l err] = waterfill_solver(P, q, C, d, max_iter)

p = diag(P);
%mu = rand(size(C,1),1);
mu = zeros(size(C,1),1);
l_aux = p.^(-1).*(C'*mu - q);
l = max(0, l_aux);

%%
    k=1;
    res = norm(C*l - d);
    while k < max_iter
        if (res > 1e-4) %this threshold is tunable, 10^-4 is good enough
            mu = mu - 0.1*(C*l - d); % gradient-descent 
            l_aux = p.^(-1).*(C'*mu - q);
            l = max(0, l_aux);
        end
        k = k+1;
        res = norm(C*l - d);
        err(k) = res;
    end
end