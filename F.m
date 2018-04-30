function [out] = F(R,T,beta,k,z)

[m,n] = size(R);
out = zeros(m*n,1);
for i = 1:m
    for j = 1:n
        new_u = BSplineTransformation([i,j],beta,k,z);
        
        out((i-1)*n+j) = BilinearApp(T,[i-new_u(1),j-new_u(2)]) - R(i,j);
    end
end
end