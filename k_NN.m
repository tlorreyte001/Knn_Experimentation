function [ kNN ] = k_NN( k, data, X, U, image, l, N, size_cls_test)
%K_NN Summary of this function goes here
%   Detailed explanation goes here
Nu=[];
for i=1:N
    Nu(i)=norm(omega(data, U, image, l, N) - omega(X, U, i, l, N));
end


Numin_indice=zeros(1,k);
[Nusort,I]=sort(Nu);
for i=1:k
    Numin_indice(i)=I(i);
end

phi=[];
for j=1:6
    a=size(intersect(Numin_indice,(size_cls_test*(j-1)+1:size_cls_test*j)));
    phi(j)=a(2);
end

[scorekNN, kNN]=max(phi);


end

