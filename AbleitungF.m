function [A,F] = AbleitungF(R,T,dxT,dyT,beta,k,z)
k1 = k(1);
k2 = k(2);
z2 = z(2);
[m,n] = size(T);
[lG,~] = size(beta);
I = zeros(n*m*32,1);
J = zeros(n*m*32,1);
V = zeros(n*m*32,1);
F = zeros(n*m,1);
Index = 0;
for i = 1:m
    t1 = floor((i-1)/k1)*z2*2;
    u = (i-1)/k1 - floor((i-1)/k1);
    
    zz(1) = ((1-u)^3)/6;
    zz(2) = (3*u^3 -6*u^2 +4)/6;
    zz(3) = (-3*u^3 + 3*u^2 +3*u+1)/6;
    zz(4) = (u^3)/6;
    for j = 1:n
        t2 = floor((j-1)/k2)*2 +1;
        u2 = (j-1)/k2 - floor((j-1)/k2);
        zz2(1) = ((1-u2)^3)/6;
        zz2(2) = (3*u2^3 -6*u2^2 +4)/6;
        zz2(3) = (-3*u2^3 + 3*u2^2 +3*u2+1)/6;
        zz2(4) = (u2^3)/6;
        
        %%BSplineTransformation
        tt   = t2 +t1;
        tt1  = t2 +t1+2;
        tt2  = t2 +t1+4;
        tt3  = t2 +t1+6;
        tt4  = t2 +t1+2*z2;
        tt5  = t2 +t1+2+2*z2;
        tt6  = t2 +t1+4+2*z2;
        tt7  = t2 +t1+6+2*z2;
        tt8  = t2 +t1+4*z2;
        tt9  = t2 +t1+2+4*z2;
        tt10 = t2 +t1+4+4*z2;
        tt11 = t2 +t1+6+4*z2;
        tt12 = t2 +t1+6*z2;
        tt13 = t2 +t1+2+6*z2;
        tt14 = t2 +t1+4+6*z2;
        tt15 = t2 +t1+6+6*z2;

        new_u = beta(tt:tt+1)*zz(1)*zz2(1);
        new_u = new_u + beta(tt1:tt1+1)*zz(1)*zz2(2);
        new_u = new_u + beta(tt2:tt2+1)*zz(1)*zz2(3);
        new_u = new_u + beta(tt3:tt3+1)*zz(1)*zz2(4);
        new_u = new_u + beta(tt4:tt4+1)*zz(2)*zz2(1);
        new_u = new_u + beta(tt5:tt5+1)*zz(2)*zz2(2);
        new_u = new_u + beta(tt6:tt6+1)*zz(2)*zz2(3);
        new_u = new_u + beta(tt7:tt7+1)*zz(2)*zz2(4);
        new_u = new_u + beta(tt8:tt8+1)*zz(3)*zz2(1);
        new_u = new_u + beta(tt9:tt9+1)*zz(3)*zz2(2);
        new_u = new_u + beta(tt10:tt10+1)*zz(3)*zz2(3);
        new_u = new_u + beta(tt11:tt11+1)*zz(3)*zz2(4);
        new_u = new_u + beta(tt12:tt12+1)*zz(4)*zz2(1);
        new_u = new_u + beta(tt13:tt13+1)*zz(4)*zz2(2);
        new_u = new_u + beta(tt14:tt14+1)*zz(4)*zz2(3);
        new_u = new_u + beta(tt15:tt15+1)*zz(4)*zz2(4);
        
        
        %%
        %new_u = BSplineTransformation([i,j],beta,k,z);
        %4 Beta's ungleich null pro Summe(jedes Beta besteht aus 2 Komponenten nach denen Abgeleitet wird)
        %Da die Ableitung von u nur jeweils in einer Kompnente ungleich
        %null ist
        ii = i-new_u(1);
        jj = j-new_u(2);
        a1 = BilinearApp(dyT,[ii;jj]);
        %Ableitung nach der zweiten Komponente
        a2 = BilinearApp(dxT,[ii;jj]);
        
        for o = 0:3
%           b(1) = BSpline(o,k(1),i);           
            for w = 0:3
%                 b(2) = BSpline(w,k(2),j);
                newa = a1*zz(o+1)*zz2(w+1);
                if abs(newa) > 10^-3
                    Index = Index +1;
                    I(Index) = (i-1)*n+j;
                    J(Index) = t1+t2 +2*z2*o+ 2*w;
                    V(Index) = -newa;%a1*zz(o+1)*zz2(w+1);
                end
%               A((i-1)*n+j, t(1)+t(2) +2*z(2)*o+ 2*w)  = -a(1)*b(1)*b(2);

                newa = a2*zz(o+1)*zz2(w+1);
                if abs(newa) > 10^-3
                    Index = Index +1;
                    I(Index) = (i-1)*n+j;
                    J(Index) = t1+t2 +2*z2*o+ 2*w +1;
                    V(Index) = -newa;%a2*zz(o+1)*zz2(w+1);
%                   A((i-1)*n+j, t(1)+t(2) +2*z(2)*o+ 2*w +1)  = -a(2)*b(1)*b(2);
                end
            end
        end

        
        F((i-1)*n+j) = BilinearApp(T,[i-new_u(1),j-new_u(2)]) - R(i,j);
    end
end
I = I(1:Index);
J = J(1:Index);
V = V(1:Index);
A = sparse(I,J,V,n*m,lG);
end