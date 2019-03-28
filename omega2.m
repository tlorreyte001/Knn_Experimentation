function [ omega2 ] = omega2( X, omeg, MUi, SIGMA )
%OMEGA2 Summary of this function goes here
%   Detailed explanation goes here

omega2=(sqrt(1/(2*pi*det(SIGMA))))*exp((-1/2)*(omeg-MUi)'*SIGMA^(-1)*(omeg-MUi);

end

