function [Xpredicted, Ppredicted, Pcrosscovariance] = UT(X_Previous, Pxx_Previous, f, param_f, Weight, Sigma_Rep, Measur_Dim)
%% Description
%Implementation of the Unscented Transformation function

%Author: Victor Araujo Vieira.
%e-mail: icevct@gmail.com
%University of Brasilia - Brazil.

%Inputs
%X_Previous: Initial State matrix
%Pxx_Previous: Covariance matrix
%f: The state function
%param_f: The state's function parameters
%Weight: The weight that was handled to the main UKF function
%Sigma_Rep: The option parameter, which chooses the sigma representation
%           that will be used
%Measur_Dim: The amount of rows in the measurements function

%Outputs
%Xpredicted: The predicted state
%Ppredicted: The predicted covariance matrix
%Pcrosscovariance: The cross covariance matrix

%% Implementation

%n: amount of states or measurements
if(nargin < 7)
    n = size(X_Previous, 1);
else
    n = Measur_Dim;
end

%Functions that computes the sigma points
%Sigma_Rep its a option variable, and it will choose the sigma
%representation that will be used
switch Sigma_Rep
    case 0
        [SigmaPointsBefore, Weights] = EvenHomiSySigRep(X_Previous, Pxx_Previous);
    case 1
        [SigmaPointsBefore, Weights] = SigRepJulier1995(X_Previous, Pxx_Previous, Weight);
    case 2
        [SigmaPointsBefore, Weights] = HoMiSySigRep(X_Previous, Pxx_Previous, Weight);
    case 3
        [SigmaPointsBefore, Weights] = RhoMiSigRep(X_Previous, Pxx_Previous, Weight);
    otherwise
         error('There are not that option for Sigma Representation!');
end

N = size(SigmaPointsBefore,2); %n amount of sigma points rows and weights

%Generates the current sample sigma pointsm evaluating the previous sigma
%points computed by the sigma representation function
SigmaPointsNow = zeros(n, N);
for i=1:N
    if(isempty(param_f))
        SigmaPointsNow(:, i) = feval(f , SigmaPointsBefore(:, i));
    else
        SigmaPointsNow(:, i) = feval(f , SigmaPointsBefore(:, i), param_f);
    end
end

%Function that computes the XPredicted and the Ppredicted
[Xpredicted, Ppredicted] = State_or_MeasurePredicted_CovMatrixPredicted(SigmaPointsNow, Weights);

%Function that computes the cross covariance matrix Pxy
Pcrosscovariance = Pcrosscovariance_withWeights(SigmaPointsBefore, SigmaPointsNow, Weights);

end