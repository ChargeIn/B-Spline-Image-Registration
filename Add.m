
function [A] = Add(A,lambda)

[m,n] = size(A);

for i = 1:m
    A(i,i) = A(i,i)+ lambda;
end

end
