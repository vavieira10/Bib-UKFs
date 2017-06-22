function [SigmaPointsBefore, Weights] = SigRepJulier1995(State, Cov_Matrix, Weight)
%% Description
%Function that computes the sigma points according to
%SigmaRepresentationJulier1995. (Symmetric set of [1])
%Paper Table I, Sigma representation 1 

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

n = size(State, 1); % amount of states
N = 2*n + 1; % amount of sigma points
Weights = zeros(N, 1); 
SigmaPointsBefore = zeros(n, N);

if Weight <= -n
    error('Weight must have a value higher than -n!\n');
end

NewWeight = Weight/(n + Weight);
% loop that computes the weights
for i = 1:N
    Weights(i, 1) = 1/(2*(n + Weight));
end
Weights(1, 1) = NewWeight; % first weight must be the initial weight passed as argument

DeltaSigmaPoint = chol((n + Weight)*Cov_Matrix,'lower');

SigmaPointsBefore(:, 1) = State; % the first sigma point is the X previous
% loop that computes the sigma points
for i = 1:n
   SigmaPointsBefore(:, i + 1) = State + DeltaSigmaPoint(:, i);
   SigmaPointsBefore(:, i + n + 1) = State - DeltaSigmaPoint(:, i);
end

end

