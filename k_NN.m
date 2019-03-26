function [ Sortie ] = k_NN( k, adr, W, l, U, N )
%K_NN Summary of this function goes here
%   Detailed explanation goes here

Sortie=[];

%Extraction des données de test
fld = dir(adr);
nb_elt_test = length(fld);
for i=1:nb_elt_test
    if fld(i).isdir == false
        img = double(imread([adr fld(i).name]));
        data_test = img(:);
        
        %Calcul d'omega pour l'image test
        omega=[];
        for j=1:l
            omega=[omega, data_test.'*U(:,N-j+1)];
        end
        
        %Calcul de Nu
        Nu=[];
        for j=1:N
            Nu=[Nu, norm(omega-W(j,:))];
        end
        
        Numin_indice=[];
        [Nusort,I]=sort(Nu);
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
        
        Sortie=[Sortie kNN];
    end
end


