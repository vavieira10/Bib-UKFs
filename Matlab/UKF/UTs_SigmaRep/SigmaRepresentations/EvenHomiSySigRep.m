function [SigmaPointsBefore, Weights] = EvenHomiSySigRep(State, Cov_Matrix)
%% Description
%Function that computes the sigma points according to
%EvenHomiSySigRep. 
%It's the Paper Table I, Sigma representation 2, but with no Weight

%Author: Victor Araujo Vieira.
%e-mail: icevct@gmail.com
%University of Brasilia - Brazil.

%Inputs
%State: Initial State matrix
%Cov_Matrix: Covariance matrix

%Outputs
%SigmaPointsBefore: Sigma Points that was generated
%Weights: Vector of weights 

%% Implementation

n = size(State, 1); % amount of states
SigmaPointsBefore = zeros(n, 2*n);

DeltaSigmaPoint = (chol(n*Cov_Matrix))';
Weights = repmat(1/(2*n), 2*n, 1); % computes the weights vector

% loop that computes the sigma points
for i = 1:n
   SigmaPointsBefore(:, i) = State + DeltaSigmaPoint(:, i);
   SigmaPointsBefore(:, i + n) = State - DeltaSigmaPoint(:, i);
end

