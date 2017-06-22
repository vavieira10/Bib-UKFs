function [SigmaPointsBefore, Weights] = RhoMiSigRep(State, Cov_Matrix, Weight)
%% Descriptions
%Function that computes the sigma points according to
%RhoMinimimSigmaRepresentation. (Minimum set of [12])
%Paper Table I, Sigma representation 6 

%Author: Victor Araujo Vieira.
%e-mail: icevct@gmail.com
%University of Brasilia - Brazil.

%Inputs
%State: Initial State matrix
%Cov_Matrix: Covariance matrix
%Weight: The weight that was handled to the main UKF function

%Outputs
%SigmaPointsBefore: Sigma Points that was generated
%Weights: Vector of weights 

%% Implementation
% most of the variable names follows the paper implementation
n = size(State, 1); % number of states
N = n + 1; % number of sigma points
Weights = zeros(N, 1); 
RhoMatrix = ones(n);
OneVector = ones(n, 1);
IdentityMatrix = eye(n);
SigmaPointsBefore = zeros(n, N);

if (Weight < 0 || Weight > 1)
    error('Weight must have a value higher than 0 or lesser than 1!');
end

Rho = sqrt((1 - Weight)/n); % the rho weight
RhoMatrix = (Rho^2)*RhoMatrix; 
C = chol(IdentityMatrix - RhoMatrix)'; % the C matrix, according to the paper
AuxiliarMatrix = ((C^-1)*Weight*RhoMatrix*((C')^-1));

for i = 1:n
    Weights(i, 1) = AuxiliarMatrix(i, i);
end
W = diag(Weights(1:n));
Weights(N, 1) = Weight;

AuxMatrix1 = (chol(Cov_Matrix)')*C*(((chol(W)'))^-1);
AuxMatrix2 = (-Rho)*chol(Cov_Matrix)'*(OneVector/sqrt(Weight));
FinalMatrix = [AuxMatrix1, AuxMatrix2];

for i = 1:N
    SigmaPointsBefore(:, i) = FinalMatrix(:, i) + State;
end

end

