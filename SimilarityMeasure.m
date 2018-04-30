function [out] = SimilarityMeasure(R,T,beta,k,z)
out = 0;
k1 = k(1);
k2 = k(2);
z2 = z(2);
[m,n] = size(R);


%outw = (norm(F(R,T,beta,k,z),2))^2;
%outw = outw/2;

for i = 1:m
    y1 = floor((i-1)/k2)*z2*2;
    u = (i-1)/k1 - floor((i-1)/k1);
    zz1 = ((1-u)^3)/6;
    zz2 = (3*u^3 -6*u^2 +4)/6;
    zz3 = (-3*u^3 + 3*u^2 +3*u+1)/6;
    zz4 = (u^3)/6;
    for j = 1:n
%       new_u = BSplineTransformation([i,j],beta,k,z);
        y2 = floor((j-1)/k2)*2 +1;
        u2 = (j-1)/k2 - floor((j-1)/k2);
        
        zz21 = ((1-u2)^3)/6;
        zz22 = (3*u2^3 -6*u2^2 +4)/6;
        zz23 = (-3*u2^3 + 3*u2^2 +3*u2+1)/6;
        zz24 = (u2^3)/6;

        t   = y2 +y1;
        t1  = y2 +y1+2;
        t2  = y2 +y1+4;
        t3  = y2 +y1+6;
        t4  = y2 +y1+2*z2;
        t5  = y2 +y1+2+2*z2;
        t6  = y2 +y1+4+2*z2;
        t7  = y2 +y1+6+2*z2;
        t8  = y2 +y1+4*z2;
        t9  = y2 +y1+2+4*z2;
        t10 = y2 +y1+4+4*z2;
        t11 = y2 +y1+6+4*z2;
        t12 = y2 +y1+6*z2;
        t13 = y2 +y1+2+6*z2;
        t14 = y2 +y1+4+6*z2;
        t15 = y2 +y1+6+6*z2;


        new_u = beta(t:t+1)*zz1*zz21;
        new_u = new_u + beta(t1:t1+1)*zz1*zz22;
        new_u = new_u + beta(t2:t2+1)*zz1*zz23;
        new_u = new_u + beta(t3:t3+1)*zz1*zz24;
        new_u = new_u + beta(t4:t4+1)*zz2*zz21;
        new_u = new_u + beta(t5:t5+1)*zz2*zz22;
        new_u = new_u + beta(t6:t6+1)*zz2*zz23;
        new_u = new_u + beta(t7:t7+1)*zz2*zz24;
        new_u = new_u + beta(t8:t8+1)*zz3*zz21;
        new_u = new_u + beta(t9:t9+1)*zz3*zz22;
        new_u = new_u + beta(t10:t10+1)*zz3*zz23;
        new_u = new_u + beta(t11:t11+1)*zz3*zz24;
        new_u = new_u + beta(t12:t12+1)*zz4*zz21;
        new_u = new_u + beta(t13:t13+1)*zz4*zz22;
        new_u = new_u + beta(t14:t14+1)*zz4*zz23;
        new_u = new_u + beta(t15:t15+1)*zz4*zz24;
        out = out + (BilinearApp(T,[i-new_u(1),j-new_u(2)])-R(i,j))^2;
    end
end
out = out/2;
end