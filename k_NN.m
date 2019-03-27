function [ classes ] = k_NN( k, adr, X, l, U, N, P )
%Arguments
% k : Hyperparamètre k du classifieur
% adr : Adresse du dossier contenant les images de test
% X : Données de l'apprentissage
% l : Dimmension du facespace
% U : Eigenfaces
% N : Dimmension des données d'apprentissage
% P : Nombre de points pour une image
%%%%%

%Sortie
% classes : Vecteur contenant les classes des images du dossier de test
%%%%%

classes=[];

%Extraction des données de test
fld = dir(adr);
nb_elt_test = length(fld);
data_test=[];
for i=1:nb_elt_test
    if fld(i).isdir == false
        img = double(imread([adr fld(i).name]));
        data_test = [data_test, img(:)];
    end
end

%Calcul de la moyenne des données
xbarre=zeros(P,1);
for i=1:nb_elt_test-2
    xbarre=xbarre+data_test(:,i);
end
xbarre=xbarre./nb_elt_test;

%Retranchage de la moyenne
for i=1:nb_elt_test-2
    data_test(:,i)=(data_test(:,i)-xbarre);
end

for i=1:nb_elt_test-2 %Boucle sur tous les visages des données de test
    
    %Calcul de Nu
    Nu=[];
    for j=1:N
        Nu(i,j)=norm(omega(data_test,U,i,l, N)-omega(X,U,j,l,N));
    end

    Numin_indice=[];
    [Nusort,I]=sort(Nu(i,:));
    for j=1:k
        Numin_indice=[Numin_indice, I(j)];
    end

    %Calcul de phi
    phi=[];
    for j=1:6
        a=size(intersect(Numin_indice,(((j-1)*10+1:j*10))));
        phi(j)=a(2);
    end
    
    [scorekNN, kNN]=max(phi);
    classes=[classes kNN];
end
    



