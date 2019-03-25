function [ omega ] = omega( X, U, personne, l, N )
%OMEGA Summary of this function goes here
%   Detailed explanation goes here
omega=[];
for i=1:l
    omega(i)=X(:,personne).'*U(:,N-i+1);
end

end

