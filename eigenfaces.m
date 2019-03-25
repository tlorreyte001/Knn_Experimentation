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
F = zeros(192*Nc,168*max(size_cls_trn));
for i=1:Nc
   for j=1:size_cls_trn(i)
         pos = sum(size_cls_trn(1:i-1))+j;
         F(192*(i-1)+1:192*i,168*(j-1)+1:168*j) = reshape(data_trn(:,pos),[192,168]);
   end
end
figure;
imagesc(F);
colormap(gray);
axis off;


%% ACP
%Question 1 : Calcul des vecteurs propres

X=data_trn;
Xmoy=sum(X,2)/N;
for i=1:N
    X(:,i)=(X(:,i)-Xmoy)/sqrt(N);
end
[V,lambda]=eig(X'*X);
U=X*V*(V'*X'*X*V)^-(1/2);

%Question 2 :
% figure;
% for i=1:N
%     F=reshape(U(:,i),[192,168]);
%     imagesc(F);
%     axis off;
%     subplot(6,10,i);
% end
% colormap(gray);



%Question 3 : 
% figure;
personne=19;
% Xpro=zeros(P,N);
% for l=1:N-1
%     subplot(6,10,l)
%     Xpro(:,l)=U(:,N).'*(X(:,personne))*U(:,N);
%     for i=2:l
%         Xpro(:,l)=Xpro(:,l)+U(:,N-i).'*(X(:,personne))*U(:,N-i);
%     end
%     %F=reshape(Xpro,[192,168]); %Sans la moyenne 
%     F=reshape(Xpro(:,l)+Xmoy/sqrt(N),[192,168]);
%     imagesc(F);
%     colormap(gray);
%     axis off;
% end

%Question 4 : 
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

%% CLASSIFICATION

%Question 1


k=7;
Nu=[];
for i=1:N
    Nu(i)=norm(omega(X,U,personne,l, N)-omega(X,U,i,l,N));
end


Numin_indice=zeros(1,k);
[Nusort,I]=sort(Nu);
for i=1:k
    Numin_indice(i)=I(i);
end

phi=[];
for j=1:6
    a=size(intersect(Numin_indice,(size_cls_trn(1)*(j-1)+1:size_cls_trn(1)*j)));
    phi(j)=a(2);
end

[scorekNN, kNN]=max(phi);

%% Matrice de confusion

% Data extraction
% Test set
adr = './database/test1/';
fld = dir(adr);
nb_elt_test = length(fld);
% Data matrix containing the testing images in its columns 
data_test = []; 
% Vector containing the class of each testing image
lb_test = []; 
for i=1:nb_elt_test
    if fld(i).isdir == false
        lb_test = [lb_test ; str2num(fld(i).name(6:7))];
        img = double(imread([adr fld(i).name]));
        data_test = [data_test img(:)];
    end
end
% Size of the testing set
[P_test,size_test] = size(data_test);
% Classes contained in the testing set
[~,I]=sort(lb_test);
data_test = data_test(:,I);
[cls_test,bd,~] = unique(lb_test);
Nc = length(cls_test); 
% Number of training images in each class
size_cls_test = [bd(2:Nc)-bd(1:Nc-1);N-bd(Nc)+1];