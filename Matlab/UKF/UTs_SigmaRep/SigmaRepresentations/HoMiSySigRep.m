function [SigmaPointsBefore, Weights] = HoMiSySigRep(State, Cov_Matrix, Weight)
%% Description
%Function that computes the sigma points according to
%HoMinimimSymmetricSigmaRepresentation. (Symmetric set of [7])
%Paper Table I, Sigma representation 2 

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

if Weight > 1
    error('Weight must have a value lesser than 1!\n');
end

n = size(State, 1); % amount of states
N = 2*n + 1; % amount of sigma points
Weights = zeros(N, 1); 
SigmaPointsBefore = zeros(n, N);

% loop that generates the weights
for i = 1:N
    Weights(i, 1) = (1 - Weight)/(2*n);
end

Weights(1, 1) = Weight; % first weight must be the initial weight passed as argument

DeltaSigmaPoint = chol((n/(1 - Weight))*Cov_Matrix,'lower');

SigmaPointsBefore(:, 1) = State; % the first sigma point is the X previous
% loop that computes the sigma points
for i = 1:n
   SigmaPointsBefore(:, i + 1) = State + DeltaSigmaPoint(:, i);
   SigmaPointsBefore(:, i + n + 1) = State - DeltaSigmaPoint(:, i);
end




