% P. Vallet (Bordeaux INP), 2019

clc;
clear all;
close all;

%% Data extraction
% Training set
adr = './database/training1/';
fld = dir(adr);
nb_elt = length(fld);
% Data matrix containing the training images in its columns 
data_trn = []; 
% Vector containing the class of each training image
lb_trn = []; 
for i=1:nb_elt
    if fld(i).isdir == false
        lb_trn = [lb_trn ; str2num(fld(i).name(6:7))];
        img = double(imread([adr fld(i).name]));
        data_trn = [data_trn img(:)];
    end
end
% Size of the training set
[P,N] = size(data_trn);
% Classes contained in the training set
[~,I]=sort(lb_trn);
data_trn = data_trn(:,I);
[cls_trn,bd,~] = unique(lb_trn);
Nc = length(cls_trn); 
% Number of training images in each class
size_cls_trn = [bd(2:Nc)-bd(1:Nc-1);N-bd(Nc)+1]; 
% Display the database
%F = zeros(192*Nc,168*max(size_cls_trn));
%for i=1:Nc
%    for j=1:size_cls_trn(i)
%          pos = sum(size_cls_trn(1:i-1))+j;
%          F(192*(i-1)+1:192*i,168*(j-1)+1:168*j) = reshape(data_trn(:,pos),[192,168]);
%    end
%end
%figure;
%imagesc(F);
%colormap(gray);
%axis off;

%% ACP
%Calcul des vecteurs propres de R
xbarre=zeros(P,1);
for i=1:N
    xbarre=xbarre+data_trn(:,i);
end
xbarre=xbarre./N;
X=data_trn;
for i=1:N
    X(:,i)=(1/sqrt(N)).*(X(:,i)-xbarre);
end
[V,lambda]=eig(X'*X);
U=X*V*(lambda)^(-0.5);
% U=X*V(V'*X'*X*V)^(-1/2);

%Représentation des n eigenfaces
% figure;
% for i=1:N
%     F=reshape(real(U(:,i)),[192,168]);
%     subplot(5,6,i);
%     imagesc(F);
%     axis off;
% end
% colormap(gray);

%Image reconstruite
% personne=1;
% figure;
% Xpro=zeros(P,N);
% for l=1:N-1
%     subplot(6,10,l);
%     Xpro(:,l)=U(:,N).'*(X(:,personne))*U(:,N);
%     for i=2:l
%         Xpro(:,l)=Xpro(:,l)+U(:,N-i).'*(X(:,personne))*U(:,N-i);
%     end
%     F=reshape(Xpro(:,l)+xbarre/sqrt(N),[192,168]);
%     imagesc(F);
%     colormap(gray);
%     axis off;
% end

%Ratio de reconstruction
% figure;
ratio=[];
ratio(1)=lambda(N,N)/sum(diag(lambda));
for l=2:N
    ratio(l)=ratio(l-1)+lambda(N+1-l,N+1-l)/sum(diag(lambda));
end
% plot(ratio);
% title('Ratio de reconstruction');
l=1;
while ratio(l)<0.90
    l=l+1;
end

%% Classification

% k=15;
% lb_test=k_NN(k, './database/test6/', X, l, U, N, P)
% lb=[ones(1,12), 2*ones(1,12) ,3*ones(1,12) ,4*ones(1,12),5*ones(1,12), 6*ones(1,12)];
% 
% [C,err_rate]=confmat(lb',lb_test');
% C
% err_rate

%% Classification Gaussienne

%Calcul des moyennes intra-classe à partir des données d'apprentissage
MU=[];
Cj=size_cls_trn(1);
nbre_classes=size(size_cls_trn);
for j=1:nbre_classes
    o=zeros(1,l);
    for i=1:Cj
        o=o+omega(X, U, i+((j-1)*10), l, N);
    end
    MU=[MU;o];
end
MU=(1/Cj).*MU;

%Calcul de la covariance
SIGMA=zeros(l,l);
for j=1:nbre_classes
    for i=1:Cj
        SIGMA=SIGMA+(omega(X, U, i+((j-1)*10), l, N))'*(omega(X, U, i+((j-1)*10), l, N));
    end
end

%Question 1



    