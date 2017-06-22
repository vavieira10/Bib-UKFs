Unscentend Kalman Filter Library - Python

Author: Victor Araujo Vieira.
Advisor: Henrique Menegaz
e-mail: victor.4.vieira@gmail.com
University of Brasilia - Brazil.
Version: 1.0
-----------------------------------------------------------------------------------------------------------------------------

-Sigma representations implemented

Sigma representation of Julier 1995 (Table I, 1), implemented as sig_rep_julier_1995;
Homogenemous Minumum Symmetric Sigma Representation(Corollary 4), implemented as ho_mi_sy_sig_rep;
Rho Minimum Sigma Representation(Corollary 5), implemented as rho_mi_sig_rep;
Even Minumum Symmetric Sigma Representation(Corollary 3), implemented as even_homi_sy_sig_rep;

-----------------------------------------------------------------------------------------------------------------------------

-Software required for using the library: Python3 (The implementation was made in Python 3.5.2), NumPy 1.12.1 (Version 									   used), SciPy 0.19.0 (Version used), Matplotlib 1.5.1 (Optional for running 										   the sample code pendulum.py), git (Optional)

-Steps for using the library:
		-Move the UKF folder to your project folder;
		-If the python file you're going to use the library is in the same directory as the UKF folder, then:
			- Add the following lines at the beggining of your code: 
				import sys 
				sys.path.append('./UKF')
				from Ad_UKF import *
		-Else, I suggest that you follow the steps bellow:
			-At the python file you're going to use the library, add the following lines at the beggining of your code:
				import sys 
				sys.path.append('pathToLibraryFolder')
				Where pathToLibraryFolder is the path to the UKF folder in your computer, based on your current folder;
		-Now, the library can be used, read the ad_ukf main function description below;

-----------------------------------------------------------------------------------------------------------------------------

Additive UKF dynamic system functions

	Xk = f(Xk-1, k) + Qk (State)
	Yk = h(Xk , k) + Rk (Measurements)


Function definition
- ad_ukf(x_previous, pxx_previous, state_funct, param_state_funct, measur_funct, params_measur_funct, Q, R,
		   measurements, sig_rep=0, *args)

Inputs
	
	-ALL OF THE VECTORS AND MATRICES REFERENCED BELOW, ARE NUMPY TYPES, THAT IS, A numpy.array VARIABLE WITH THE SPECIFIED DIMENSIONS

	-x_previous: 
		-Object type: A column vector, Nx1;
		-Description: Initial state estimate;

	-pxx_previous: 
		-Object type: A square matrix, NxN;
		-Description: Initial covariance matrix;

	-state_funct: 
		-Object type: A python function;
		-Description: The definiton of the state process function, which is a python function;

	-param_state_funct: 
		-Object type: A row vector, 1xA, where A is the amount of arguments;
		-Description: The state's function parameters;

	-measur_funct: 
		-Object type: A python function;
		-Description: The definiton of the measurements function, which is a python function;

	-params_measur_funct: 
		-Object type: A row vector, 1xA, where A is the amount of arguments;
		-Description: The measurements function parameters;

	-Q: 
		-Object type: A scalar, vector or matrix;
		-Description: State noise matrix;

	-R: 
		-Object type: A scalar, vector or matrix;
		-Description: Measurements noise matrix;

	-measurements:
		-Object type: A column or row vector, Mx1 or 1xM, where M is the amount of measurements;
		-Description: The measurments samples;

	-sig_rep: 
		-Object type: A scalar, in the range [0, 3]. 
		-Description: The option parameter, which chooses the sigma representation that will be used. If no option is 				  handled to the function, the default sigma representation will be EvenHomiSySigRep(sig_rep = 0);

	-*args (Sigma Representation Tunning Parameter): 
		-Object type: A scalar;
		-Description: The weight that is used in the Sigma Representation functions. No weight is used in the 						 even_homi_sy_sig_rep sigma representation, so there's no need to handle any argument.

Outputs

	-x_posterior: 
		-Object type: A column vector, Nx1;
		-Description: Corrected state matrix;

	-p_posterior: 
		-Object type: A square matrix, NxN;
		-Description: Corrected covariance matrix;

-----------------------------------------------------------------------------------------------------------------------------

An example of how to use:

import numpy as np
import math as mth

DT = 0.01
Q = 0.01*np.array([[(DT**3)/3, (DT**2)/2], [(DT**2)/2, DT]])
R = 0.1
param_F = np.zeros((2))
param_F[0] = DT
param_F[1] = 9.81
y[0, 0] = mth.sin(x[0, 0]) + mth.sqrt(R)*np.random.randn(1, 1) # measure

x_previous = np.array([[1.6], [0]]) # 2x1 vector
pxx_previous = 0.1*np.eye(2)

x_posterior, p_posterior = ad_ukf(x_previous, pxx_previous, state_function, param_F, measurement_function, [], Q, R,
					  			  y[0, 0, 2, 0.2)


#State function
def state_function(x, params):
	res = np.zeros((2, 1))

	res[0, 0] = x[0, 0] + x[1, 0]*params[0]
	res[1, 0] = x[1, 0] - params[1]*mth.sin(x[0, 0])*params[0]

	return res

#Measurement function
def measurement_function(x):
	res = np.zeros((1, 1))

	res[0, 0] = mth.sin(x[0, 0])

	return res

-----------------------------------------------------------------------------------------------------------------------------