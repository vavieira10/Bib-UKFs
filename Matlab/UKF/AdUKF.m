function [X_Posterior, Pxx_Posterior] = AdUKF(X_Previous, Pxx_Previous, f, param_f, h, param_h, Q, R, Measurements, Sigma_Rep, Sigma_Rep_Tuning_Param)
%% Description
%The main Addictive UKF function. It calls the UT and the UKF correction
%equation

%Author: Victor Araujo Vieira.
%e-mail: icevct@gmail.com
%University of Brasilia - Brazil.

%Inputs
%X_Previous: Initial State matrix
%Pxx_Previous: Covariance matrix
%f: The state function
%param_f: The state's function parameters
%h: The measurements function
%param_h: The measurements function parameters
%Q: State noise matrix 
%R: Measurements noise matrix
%Measurements: The measurments samples
%Sigma_Rep: The option parameter, which chooses the sigma representation
%           that will be used
%Sigma_Rep_Tuning_Param: The weight that was handled to the main UKF function

%Outputs
%X_Posterior: Corrected state matrix
%Pxx_Posterior: Corrected covariance matrix

%% Handling default values and initializing variables
%defining the default values for Sigma_Rep and Weight
if(isempty(Sigma_Rep))
    Sigma_Rep = 0; % the default sigma representation its EvenHomiSySigRep
    Sigma_Rep_Tuning_Param = [];
end

%Assuming that there's a better default value for every sigma
%representation function, the following piece of code will define the best
%for each one
if(nargin < 11 && Sigma_Rep ~= 0)
     error('No Sigma Representation Turning Parameter detected!');
end
y_dimension = size(Measurements, 1); % amount of measurements

%% State
[Xpredicted, Ppredicted] = UT(X_Previous, Pxx_Previous, f, param_f, Sigma_Rep_Tuning_Param, Sigma_Rep);
Ppredicted = Ppredicted + Q;

%% Measurement 
[Ypredicted, Pmeasurement, Pcrosscovariance] = UT(Xpredicted, Ppredicted, h, param_h, Sigma_Rep_Tuning_Param, Sigma_Rep, y_dimension);
Pmeasurement = Pmeasurement + R;

%% Kalman Filter correction
[X_Posterior, Pxx_Posterior] = KF_Correction(Xpredicted, Ppredicted, Ypredicted, Pmeasurement, Pcrosscovariance, Measurements);

