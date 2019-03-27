function [ omega ] = omega( X, U, image, l, N )
%Arguments
% X : Matrice de données à traiter
% U : Eigenfaces
% image : indice de l'image à traiter
% l : Dimmension du facespace
% N : Dimmension des données d'apprentissage
%%%%%

%Sortie
% omega : Vecteur des l composantes principales
%%%%%

omega=[];
for i=1:l
    omega(i)=X(:,image).'*U(:,N-i+1);
end

end

