function [beta,regi,Dssd] =  ImagePyramid(R,T,k)
%R = vorlage , T = verschobenes Bild, k = [k1,k2,...] Schrittweite 
%Hier 2D-Model daher k = [k1,k2]
%test ob Dimension von k stimmt
%  beta = [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5]';
%  beta = expand(beta,[7,7],[5,5])
%  length(beta)
test = length(k);
if test ~= 2 
    error('k hat nicht genug oder zuviele Argumente');
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


test = max(n(1),n(2));
%Berechnen der groesse der Pyramide des Bildes
ebene = 1;
while(test > 64)
    test = test/2;
    ebene =ebene+1;
end
%Berechnen der groesse der Parameter-Pyramide
test = max(floor(n(1)/2),floor(n(2)/2));
P_ebene = 1;
Max = max(n(1),n(2));
Max_k = max(k(1),k(2));
while(test > Max_k);
    P_ebene = P_ebene +1;
    test = floor(test/2);
end
P_ebene
P_ebene_b = true;
z = [0;0];
% 5 x 5 starting mesh *2 fuer x,y Koordinate 
beta = zeros(5*5*2,1);
while(ebene > 0)
    ebene
    new_R = R;
    new_T = T;
    for i = 1:ebene-1
        new_R = impyramid(new_R,'reduce');
        new_T = impyramid(new_T,'reduce');
    end
    [n(1),n(2)] = size(new_R);
    if P_ebene_b == true
        for i = 1:P_ebene
            i
            k(1) = (n(1)/2^i);
            k(2) = (n(2)/2^i);
            z_old = z;
            for j = 1:2
                z(j) = ceil(n(j)/k(j))+3;
            end
            beta = expand(beta,z,z_old);
            [beta,regi,Dssd] = ImageRegistration(new_R,new_T,k,beta);
        end
        P_ebene_b = false;
    else
        k(1) = ceil(n(1)/(z(1)-3));
        k(2) = ceil(n(2)/(z(2)-3));
        [beta,regi,Dssd] = ImageRegistration(new_R,new_T,k,beta);
    end
    ebene = ebene-1;
end

% Image = zeros(n(1),n(2));
% for i = 1:n(1);
%     for j = 1:n(2)
%         new_u = BSplineTransformation([i,j],beta,k,z);
%         Image(i,j)=BilinearApp(T,[i-new_u(1),j-new_u(2)]);
%     end
% end
end


function [new_beta] = expand(beta,z_new,z_old)
new_beta = zeros(z_new(1)*z_new(2)*2,1);
[m,~] = size(beta);
z = 1;
t = 1;
l = false;
z2 = (z_new(2)-3)/2+3;
for i = 1:z_new(1)
    for j= 1:z2
        %Achtung die Aeusseren 3 Kontrolpunkte müssen weiterhin aussen
        %bleibn
        if i > 3 && j > 3 
            %false = alle ungeraden Zeilen (beginnend bei 1)
            if l == false;
                new_beta(z:z+1)  = beta(t:t+1);
                z = z+2;
                if mod(t+1,2*z_old(2)) ~= 0
                    new_beta(z:z+1)= (beta(t:t+1)+beta(t+2:t+3))/2; 
                else
                    new_beta(z:z+1)= beta(t:t+1);                    
                end
                z = z+2;
                if j == z2
                    t =  t - 2*z_old(2);
                end
            else
                if m-t > 2*z_old(2)
                    t2 = t+2*z_old(2);
                    new_beta(z:z+1) = (beta(t:t+1) + beta(t2:t2+1))/2;
                    z = z+2;
                    if mod(t+1,2*z_old(2)) ~= 0
                        t2 = t+2*z_old(2);
                        new_beta(z:z+1)= (beta(t:t+1)+beta(t+2:t+3)...
                                        +beta(t2:t2+1)+beta(t2+2:t2+3))/4;
                    else
                        t2 = t+2*z_old(2);
                        new_beta(z:z+1)= (beta(t:t+1)+beta(t2:t2+1))/2;
                    end
                    z = z+2;
                else
                    new_beta(z:z+1)= beta(t:t+1);
                    z = z+2;
                    if mod(t+1,2*z_old(2)) ~= 0
                        new_beta(z:z+1)= (beta(t:t+1)+beta(t+2:t+3))/2;
                    else
                        new_beta(z:z+1)= beta(t:t+1);
                    end
                    z = z+2;
                end
            end
            if j == z2
                l = ~l;
            end
        elseif i < 4 && j > 3
            new_beta(z:z+1)  = beta(t:t+1);
            z = z+2;
            if mod(t+1,2*z_old(2)) ~= 0
                new_beta(z:z+1)= (beta(t:t+1)+beta(t+2:t+3))/2; 
            else
                new_beta(z:z+1)= beta(t:t+1);                    
            end
            z = z+2;
        elseif j < 4 && i > 3
            if l == false;
                new_beta(z:z+1)  = beta(t:t+1);
                z = z+2;
            else
                if m-t > 2*z_old(2)
                    t2 = t+2*z_old(2);
                    new_beta(z:z+1) = (beta(t:t+1) + beta(t2:t2+1))/2;
                    z = z+2;
                else
                    new_beta(z:z+1)= beta(t:t+1);
                    z = z+2;
                end
            end
        else
            new_beta(z:z+1) = beta(t:t+1); 
            z = z+2;
        end
        t = t+2;
    end
end

end





