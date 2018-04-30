function [out] = DSSD(R,T,beta,k,z)
out = 0;
[m,n] = size(R);
%outw = (norm(F(R,T,beta,k,z),2))^2;
%outw = outw/2;
for i = 1:m
    for j = 1:n
        new_u = BSplineTransformation([i,j],beta,k,z);
        out = out + (BilinearApp(T,[i-new_u(1),j-new_u(2)])-R(i,j))^2;
    end
end
out = out/(m*n);

end