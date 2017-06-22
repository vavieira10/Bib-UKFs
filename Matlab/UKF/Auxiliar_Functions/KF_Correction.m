function [Xposterior, Pposterior] = KF_Correction(Xpredicted, Ppredicted, Ypredicted, Pmeasurement, Pcrosscovariance, Measurements)
%% Description
%Function that computes the KF equation

%Author: Victor Araujo Vieira.
%e-mail: icevct@gmail.com
%University of Brasilia - Brazil.

%Inputs
%Xpredicted: State predicted matrix
%Ppredicted: Predicted state covariance matrix
%Ypredicted: Measurement predicted matrix
%Pmeasurement: Predicted measurement covariance matrix
%Pcrosscovariance: Cross covariance matrix
%Measurements: The measurements sample

%Outputs
%Xposterior: Corrected state matrix
%Pposterior: Corrected covariance matrix

%% Implementation

Kk = Pcrosscovariance/Pmeasurement; 
Xposterior = Xpredicted + Kk*(Measurements - Ypredicted);
Pposterior = Ppredicted - Kk*Pmeasurement*(Kk');

end

