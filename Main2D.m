%Small image registration test (rotation & scalling)
k = [100,100];
original = imread('Knie.jpg');
original = im2double(original(1:600,1:600));
[x,y] = size(original);
original = [original,zeros(x,100)];
%original = im2double(original(101:300,101:295));
%original = ones(40,35);
%original2 = zeros(x+2*k(1),y+2*k(2));
%original2(k(1)+1:x+k(1),1+k(2):y+k(2)) = original;
%original = original2;
%clear original2
%Setting up the images
%original = zeros(60,55);
%for i = 0:4
%    for j = 1:10
%        for kl = 1:60
%            original(kl,j+10*i) = (i+1)/5;
%        end
%    end
%end
%scale
%scale = 0.71; % scaling factor
%rotated = imresize(original, scale);

% rotated
%alpha = 40; %Degree

%rotated = imrotate(rotated, alpha);
%[s1,s2] =size(rotated);
%rotated = rotated(3:s1,3:s2);


%shift
original = imtranslate(original,[100, 0]);
rotated =imtranslate(original,[-5, 0]);

% for i = 1:20+5
%    for j = 1:x
%       original(j,y-20+i) = 0;
%       if i < 6
%         rotated(j,y+i) = 0;
%       end
%    end
% end
[x,y] = size(original);

[beta,regi,Dssd] = ImageRegistration(original,rotated,k);

%Erzeugen des registrierten Bildes
z = [0;0];
z(1) = ceil(x/k(1))+3;
z(2) = ceil(y/k(2))+3;
U = zeros(x,y);
V = zeros(x,y);
for i = 1:x
    for j = 1:y
        new_u = BSplineTransformation([i,j],beta,k,z);
        U(i,j) = new_u(1);
        V(i,j) = new_u(2);
    end
end
figure
streamslice(1:y,1:x,V,U)
% rotated = insertMarker(rotated,[x/2-60,y/2-100],'plus','color','red','size',10);
% original = insertMarker(original,[x/2-60,y/2-100],'plus','color','red','size',10);

%Plotting
hold all
figure('name', 'Bildregistrierung', 'NumberTitle', 'Off');
%Original
subplot(1,2,1);
imshow(original);

%Rotated Image
subplot(1,2,2);
imshow(rotated);
figure
%After Imageregistration
imshow(full(regi));
title('Registriertes Bild');
%Differenz vorher
figure
subplot(1,2,1)
imshow(imabsdiff(original,rotated));
%Differenz nachher
subplot(1,2,2)
imshow(imabsdiff(original,full(regi)));
% 
% figure
% plot(1:size(Dssd_5,2),Dssd_5,1:size(Dssd_15,2),Dssd_15,1:size(Dssd_25,2),Dssd_25)
% xlabel('Iterationen')
% ylabel('Ähnlichkeitsmaß');
% title('Ähnlichkeitsmaß');
% %axis([0 100 0 0.06])
% legend('5p','15p','25p');