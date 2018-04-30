function [out] =  BilinearAppMod(T,x,Ableitung)
%modifizierung fuer die Ableitungsmatrix, da diese 1te- und 2te- Ableitung
%untereinander speichert, d.h. der Punkt x,y zweite Ableitung befindet sich an Stelle
%2*x und erste Ableitung an 2*x-1

% a11---- a12
%   |           |
%   |           |
% a21---- a22

[m,n] = size(T);
y11 = floor(x(1));
y12 = floor(x(2));
y21 = ceil(x(1));
y22 = ceil(x(2));
oldValue1 = y11;
oldValue2 = y21;

if Ableitung == 1
    y11 = 2*y11-1;
    y21 = 2*y21-1;
else
    y11 = 2*y11;
    y21 = 2*y21;
end

if x(1)>= 1 && x(1)<= m/2 && x(2)>= 1 && x(2)<= n
    a11 = T(y11,y12); 
    a12 = T(y11,y22);
    a21 = T(y21,y12);
    a22 = T(y21,y22);
else
    if y11 < 1 || y11 > m
        a11 = 0;
        a12 = 0;
    else
       if y12 < 1 || y12 > n
           a11 = 0;
       else
           a11 = T(y11,y12);
       end 
       if y22 < 1 || y22 > n
           a12 = 0;
       else
           a12 = T(y11,y22);
       end
    end
    
    if y21 < 1 || y21 > m
        a21 = 0;
        a22 = 0;
    else
        if y12 < 1 || y12 > n
            a21 = 0;
        else
            a21 = T(y21,y12);
        end
        if y22 < 1 || y22 > n
              a22 = 0;
        else
           a22 = T(y21,y22);
        end
    end
end


%da die schritweite immer nur eins ist musst man nicht den diverenzen
%quotieneten betrachten sondern nur den Abstand
if y21 ~= y11
    z1 = (oldValue2-x(1))*a11 + (x(1)-oldValue1)*a21;
    z2 = (oldValue2-x(1))*a12 + (x(1)-oldValue1)*a22;
else
    z1 = a21;
    z2 = a22;
end

if y22 ~= y12
    out = (y22-x(2))*z1 + (x(2)-y12)*z2;
else
    out = z2;
end

end