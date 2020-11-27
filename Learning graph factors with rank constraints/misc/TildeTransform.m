function Atilde = TildeTransform(A,m1,n1,m2,n2)
% function Atilde = TildeTransform(A,m1,n1,m2,n2)
% A is regarded as an m1xn1 block matrix with m2xn2 blocks.
% Atilde is the (m1*n1)x(m2*n2) matrix R(A).
% Adpated from GVL Page 713.

Atilde = [];
for j=1:n1
    for i=1:m1;
        Aij = A(1+(i-1)*m2:i*m2,1+(j-1)*n2:j*n2);
        Atilde = [Atilde; vec(Aij)'];
    end
end