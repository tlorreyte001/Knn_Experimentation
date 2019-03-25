function [ Xpro ] = Xpro( X,U,personne,P,N)
%XPRO Summary of this function goes here
%   Detailed explanation goes here
Xpro=zeros(P,N);
for l=1:N-1
    Xpro(:,l)=U(:,N).'*(X(:,personne))*U(:,N);
    for i=2:l
        Xpro(:,l)=Xpro(:,l)+U(:,N-i).'*(X(:,personne))*U(:,N-i);
    end
end

