function[beta,Image,Dssd] = ImageRegistration(R,T,k)
%R = vorlage , T = verschobenes Bild, k = [k1,k2,...] Schrittweite 
%Hier 2D-Model daher k = [k1,k2]

%test ob Dimension von k stimmt
test = length(k);
if test ~= 2 
    error('k hat nicht genug Argumente');
end
clear test

%n = [n1,n2] Anzahl an Pixel in x-,y-Richtung
[n(1),n(2)] = size(R);

%test ob R und T selben dimensionen haben
[n2(1),n2(2)] = size(T);
if n(1) ~= n2(1) || n(2) ~= n2(2)
    error('R und T haben nicht die selbe Dimension');
end
clear n2

%z = [z1,z2] Anzahl an Kontrollpunkten in x-,y- Richtung 
%Da das Bild von Kontrollpunkten umgeben wird, welche ausserhalb liegen +3
z = [0;0];
for i = 1:2
    z(i) = ceil(n(i)/k(i))+3;
end

%lG Anzahl an Kontrollpunkten
lG = z(1)*z(2);

%beta = Vektor aller Positionen der Kontrolpunkte
%fuer jeden Kontrollpunkt wird x-,y-Koordinate untereinander gespeichert 
%Startvektor = 0
beta = zeros(2*lG,1);
%spaeter zur while Schleife machen
lambda = 2;%16;
p = 0.3;
l = 0.7;
test = 99999;
new = 0;
Iteration = 1;
Dssd = zeros(1,100);
Dssd(Iteration) =  DSSD(R,T,beta,k,z);
[dxT, dyT] = imgradientxy(T,'central');
next = SimilarityMeasure(R,T,beta,k,z);
Steps=[1,0.5,0];
v = 2;
for Iteration = 1:15%while(2/(n(1)*n(2))*(test - new) > 0.000005)%
    Iteration
    [J,f] = AbleitungF(R,T,dxT,dyT,beta,k,z);

    JJ = J'*J;
    JJ = Add(JJ,lambda);
    s = (JJ)\((-J)'*f);
    alpha = Steps(1);
    step = 1;
    test = next;
    next = SimilarityMeasure(R,T,beta+alpha*s,k,z);
    %new = SimilarityMeasure(R,T,beta+0.1*s,k,z);
    while ( next > 0.995*test && alpha > 0)%step > 10^-5)
        step = step+1;
        alpha = Steps(step);
        %wegen rundungsfehler
        if alpha < 0.01
            next = test;
        else
            next = SimilarityMeasure(R,T,beta+alpha*s,k,z);
        end
        
    end
    next
    %r = (test - next)/(test-(1/2)*norm(f+J*s,2)^2)
    if test == next          
        r = 0
    else
        r = (test - next)/(-s'*alpha*J'*f)
    end
    
    beta = beta +alpha*s;
    if r < p %verkleiner lambda
       if test -new < 0  
            lambda = lambda*10; 
       else
            lambda = lambda/10;
       end
    else
        if r >l %vergroessere lambda
            lambda = lambda/10;
        end
    end
    Iteration = Iteration +1;
    Dssd(Iteration) =  DSSD(R,T,beta,k,z);
end
Image = zeros(n(1),n(2));

for i = 1:n(1);
    for j = 1:n(2)
        new_u = BSplineTransformation([i,j],beta,k,z);
        Image(i,j)=BilinearApp(T,[i-new_u(1),j-new_u(2)]);
    end
end
Dssd = Dssd(1:Iteration-1);
end

function [A] = Add(A,lambda)

[m,n] = size(A);

for i = 1:m
    A(i,i) = A(i,i)+ lambda;
end

end








