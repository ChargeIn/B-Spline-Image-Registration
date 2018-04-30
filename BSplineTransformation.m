function [out] = BSplineTransformation(x,beta,k,z)
%x,g,k zweidimensional, beta 2*lG x 1
%k = Gitterweite
%x = Koordiante des Pixels
%z = Anzahl der Kontrollpunkte in x-,y-Richtung
k1 = k(1);
k2 = k(2);
z2 = z(2);
% out = [0;0];

%Koordinaten umrechung fuer beta
%y = zeros(2,1);

%Da es negative Kontrollpunkte gibt, ist um +3 geschifted 
%beta speichert die Kontrollpunkte zeilenweise, d.h. die j-te Zeile
%faengt beim Indize (j-1)*z(2)*2+1 an.
y1 = floor((x(1)-1)/k(1))*z(2)*2;
%Da beta x-,y-Koordinate untereinander Speichert entspricht x2;
%+1 da Matlab bei 1 anfangt zu zaehlen
y2 = floor((x(2)-1)/k(2))*2 +1;
% [s,~] = size(beta);


%zz = BSpline(i,k(1),x(1));

%Berechnet den vereinfachten kubischen Spline
%i = i-tes Intervall
%k = Gitterweite
%x = Wert

%Reinfolge ist umgekehrt
u = (x(1)-1)/k1 - floor((x(1)-1)/k1);
u2 = (x(2)-1)/k2 - floor((x(2)-1)/k2);

zz1 = ((1-u)^3)/6;
zz2 = (3*u^3 -6*u^2 +4)/6;
zz3 = (-3*u^3 + 3*u^2 +3*u+1)/6;
zz4 = (u^3)/6;

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


out = beta(t:t+1)*zz1*zz21;
out = out + beta(t1:t1+1)*zz1*zz22;
out = out + beta(t2:t2+1)*zz1*zz23;
out = out + beta(t3:t3+1)*zz1*zz24;
out = out + beta(t4:t4+1)*zz2*zz21;
out = out + beta(t5:t5+1)*zz2*zz22;
out = out + beta(t6:t6+1)*zz2*zz23;
out = out + beta(t7:t7+1)*zz2*zz24;
out = out + beta(t8:t8+1)*zz3*zz21;
out = out + beta(t9:t9+1)*zz3*zz22;
out = out + beta(t10:t10+1)*zz3*zz23;
out = out + beta(t11:t11+1)*zz3*zz24;
out = out + beta(t12:t12+1)*zz4*zz21;
out = out + beta(t13:t13+1)*zz4*zz22;
out = out + beta(t14:t14+1)*zz4*zz23;
out = out + beta(t15:t15+1)*zz4*zz24;

% for i = 0:3              
%     for j = 0:3
%         t = y(2)+i*z(2)*2  +y(1)+j*2;
%         out = out + beta(t:t+1)*zz(i+1)*zz2(j+1);
%     end
% end

end
