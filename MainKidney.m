

k = [10,10];
Data = load('dyn_kidney.mat');
original = Data.F(:,:,1);
[n(1),n(2)] = size(original);
for i = 1:2
    z(i) = ceil(n(i)/k(i))+3;
end
Images = zeros(n(1),n(2));
%87 Bilder
zz = 0;
for i = 40:40
zz = zz+1;
shifted = Data.F(:,:,i);
% shifted = imrotate(shifted,30);
% shifted = shifted(51:350,51:350);
[beta,Images(:,:,zz)] = ImageRegistration(original,shifted,k);
end


U = zeros(n(1),n(2));
V = zeros(n(1),n(2));
for i = 1:1:n(1);
    for j = 1:1:n(2)
        new_u = BSplineTransformation([i,j],beta,k,z);
        U(i,j) = new_u(1);
        V(i,j) = new_u(2);
    end
end

figure
streamslice(1:n(1),1:n(2),V,U)

z = 0;
for i = 40:40
    figure
    subplot(3,3,1)
    shifted = Data.F(:,:,i);
%     shifted = imrotate(shifted,30);
%     shifted = shifted(51:350,51:350);
    imshow(shifted-original);
    z = z+1;
    subplot(3,3,2)
    imshow(Images(:,:,z)-original)
    subplot(3,3,4)
    imshow(shifted)
    subplot(3,3,5)
    imshow(Images(:,:,z))
    subplot(3,3,6)
    imshow(original)
end