function [SamplePredicted, Ppredicted] = State_or_MeasurePredicted_CovMatrixPredicted(SigmaPoints, Weights)
%% Descriptions
%Function that generates the posterior estimated state's matrix and the
%estimated state's covariance matrix

%Author: Victor Araujo Vieira.
%e-mail: icevct@gmail.com
%University of Brasilia - Brazil.

%Inputs
%SigmaPoints: SigmaPoints from the UT function
%Weights: The weight that was handled to the main UKF function

%Outputs
%SamplePredicted: posterior estimated state's matrix
%Ppredicted: posterior estimated state's covariance matrix.

N = size(Weights, 1); % amount of sigma points and weights

%State's or measures prediction
SamplePredicted = SigmaPoints*Weights;

%Covariance matrix prediction
sizes = size(SigmaPoints, 1); % number of sigma points lines
Ppredicted = zeros(sizes, sizes);

for i = 1:N
    Ppredicted = Ppredicted + Weights(i, 1)*(SigmaPoints(:, i) - SamplePredicted)*(SigmaPoints(:, i) - SamplePredicted)';
end

end

