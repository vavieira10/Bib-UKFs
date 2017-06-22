function SQRT_P=chol_or_SVD_inversion(P)
%% Description
%%This function is designed to obtain the square-root matrix of a matrix
%that can be either positive definite or semi-positive definite

%Author: Henrique Marra Menegaz. Modified from Tim Baily's UT function
%e-mail: henrique.menegaz@gmail.com
%Automation and Robotics Laboratory (LARA)
%University of Brasília - Brazil.

%creation:2012, febrary 09
%last update: 2012, febrary 09

%Inputs:
%       P:   positive definite or semi-positive definite matrix
%Outputs:
%       SQRT_P: square-root matrix of P
%% function
if det(P)==0
    [U,S,V]=svd(P);
    S=sqrt(S);
    SQRT_P = U*S;    
else    
    SQRT_P = chol(P)';
end