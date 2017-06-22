Unscentend Kalman Filter Library - MATLAB

Author: Victor Araujo Vieira.
Advisor: Henrique Menegaz
e-mail: victor.4.vieira@gmail.com
University of Brasilia - Brazil.
Version: 1.0
-----------------------------------------------------------------------------------------------------------------------------

-Sigma representations implemented

Sigma representation of Julier 1995 (Table I, 1), implemented as SigRepJulier1995;
Homogenemous Minumum Symmetric Sigma Representation(Corollary 4), implemented as HoMiSySigRep;
Rho Minimum Sigma Representation(Corollary 5), implemented as RhoMiSigRep;
Even Minumum Symmetric Sigma Representation(Corollary 3), implemented as EvenHomiSySigRep;

------------------------------------------------------------------------------------------------------------------------------

-Software required for using the library: MATLAB (The implementation was made in Matlab 2015a) and git (Optional)
	
	-Steps for using the library:
		-Open MATLAB;
		-In the directory path navegation (left side of the screen), navigate to the directory that you
		cloned the library;
		-Right click the library folder on the MATLAB path navegation, and choose the following options:
		Add to Path -> Selected Folder and Subfolders;
		-Now the library can be used, read the AdUKF main function description below;

-----------------------------------------------------------------------------------------------------------------------------

Additive UKF dynamic system functions

	Xk = f(Xk-1, k) + Qk (State)
	Yk = h(Xk , k) + Rk (Measurements)

Function definition
-[X_Posterior, Pxx_Posterior] = AdUKF(X_Previous, Pxx_Previous, f, param_f, h, param_h, Q, R, Measurements, Sigma_Rep, Sigma_Rep_Tuning_Param)

Inputs

	-X_Previous: 
		-Object type: A column vector, Nx1;
		-Description: Initial state estimate;

	-Pxx_Previous: 
		-Object type: A square matrix, NxN;
		-Description: Initial covariance matrix;

	-f: 
		-Object type: A matlab function;
		-Description: The definiton of the state process function, which is handled at MATLAB;

	-param_f: 
		-Object type: A column or row vector, Ax1 or 1xA, where A is the amount of arguments;
		-Description: The state's function parameters;

	-h: 
		-Object type: A matlab function;
		-Description: The definiton of the measurements function, which is handled at MATLAB;

	-param_h: 
		-Object type: A column or row vector, Ax1 or 1xA, where A is the amount of arguments;
		-Description: The measurements function parameters;

	-Q: 
		-Object type: A scalar, vector or matrix;
		-Description: State noise matrix;

	-R: 
		-Object type: A scalar, vector or matrix;
		-Description: Measurements noise matrix;

	-Measurements:
		-Object type: A column or row vector, Mx1 or 1xM, where M is the amount of measurements;
		-Description: The measurments samples;

	-Sigma_Rep: 
		-Object type: A scalar, in the range [0, 3]
		-Description: The option parameter, which chooses the sigma representation that will be used. If no option is 				  handled to the function, the default sigma representation will be EvenHomiSySigRep.

	-Sigma_Rep_Tuning_Param: 
		-Object type: A scalar;
		-Description: The weight that is used in the Sigma Representation functions. No weight is used in the 						 EvenHomiSySigRep sigma representation.

Outputs

	-X_Posterior: 
		-Object type: A column vector, Nx1;
		-Description: Corrected state matrix;

	-Pxx_Posterior: 
		-Object type: A square matrix, NxN;
		-Description: Corrected covariance matrix;

-----------------------------------------------------------------------------------------------------------------------------

An example of how to use:

DT = 0.01
Q  = 0.01*[DT^3/3 DT^2/2; DT^2/2 DT];
R  = 0.1;
param_F = zeros(2, 1); %parameters of the state function
param_F(1, 1) = DT;
param_F(2, 1) = 9.81;
y = sin(x(1)) + sqrt(R)*randn; % measure

f = @StateFunction;
h = @MeasurementFunction;

X_Previous = [1.6; 0];
Pxx_Previous = 0.1*eye(2);

[X_Posterior, Pxx_Posterior] = AdUKF(X_Previous, Pxx_Previous, f, param_F, h, [], Q, R, y, 1, 0.1);


%State function
function f = StateFunction(x, parameters)

f = zeros(2, 1);
f(1, 1) = x(1) + x(2)*parameters(1, 1);
f(2, 1) = x(2) - parameters(2, 1)*sin(x(1))*parameters(1, 1);

end

%Measurement function
function h = MeasurementFunction(x)

h = sin(x(1));

end

-----------------------------------------------------------------------------------------------------------------------------

