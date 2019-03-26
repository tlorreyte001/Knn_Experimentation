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
% F = zeros(192*Nc,168*max(size_cls_trn));
% for i=1:Nc
%    for j=1:size_cls_trn(i)
%          pos = sum(size_cls_trn(1:i-1))+j;
%          F(192*(i-1)+1:192*i,168*(j-1)+1:168*j) = reshape(data_trn(:,pos),[192,168]);
%    end
% end
% figure;
% imagesc(F);
% colormap(gray);
% axis off;

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
% personne=19;
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
%Omega des img de training
W=[];
for i=1:N
    wtmp=[];
    for j=1:l
        wtmp=[wtmp, X(:,i)'*U(:,N-j+1)];
    end
    W=[W;wtmp];
end

k=6;
% lb_test=k_NN(k, './database/test1/', W, l, U, N)
% module_Cj=7
% lb=[ones(1,7), 2*ones(1,7) ,3*ones(1,7) ,4*ones(1,7),5*ones(1,7), 6*ones(1,7)];
% [C,err_rate]=confmat(lb',lb_test');
% C

%% Classification Gaussienne

%Calcul des moyennes intra-classe à partir des données d'apprentissage
MU=[];
for i=1:6
    MU=[MU; 1/10*sum(W((i-1)*10+1:i*10,:))];
    
end

%Calcul de la covariance
SIGMA=0;
for i=1:6
    for j=1:10
        a= W(j+(10*(i-1)),:)-MU(i,:);
        SIGMA=SIGMA+(a*a');
    end    
end
SIGMA=(1/N)*SIGMA; %Bizarre

%Question 1: %TODO
omega=[];
for i=1:30
    fiy=sqrt(1/2*pi*det(SIGMA))*exp(-1/2*(data_trn(
    omega[
    


       
