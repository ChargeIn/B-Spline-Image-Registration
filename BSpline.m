function [out] = BSpline(i,k,x)
%Berechnet den vereinfachten kubischen Spline
%i = i-tes Intervall
%k = Gitterweite
%x = Wert

%Reinfolge ist umgekehrt
u = (x-1)/k - floor((x-1)/k);
if i == 3
    out = (u^3)/6;
elseif i == 2 
    out = (-3*u^3 + 3*u^2 +3*u+1)/6;
elseif i == 1
    out = (3*u^3 -6*u^2 +4)/6;
elseif i == 0
    out = ((1-u)^3)/6;
else
    out = 0;
end

end