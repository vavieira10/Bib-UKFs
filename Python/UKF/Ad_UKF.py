import numpy as np
import math as mth
from Sigma_Representations.sigma_representations import *
from UTs.UT import *
from Auxiliar_Functions.auxiliar_functions import *

#UKF main function implementation: ad_ukf
#Author: Victor Araujo Vieira.
#e-mail: icevct@gmail.com
#University of Brasilia - Brazil.

#ad_ukf function definition
#As there's so many variables, the arrays dimensions were not specified
#Inputs: x_previous (Array)
#		 pxx_previous (Array)
# 		 state_funct (A function. The state or measurement state function)
#		 param_state_funct (An array or an empty list [])
#		 measur_funct (A function. The state or measurement state function)
#		 params_measur_funct (An array or an empty list [])
#		 Q (Array)
# 		 R (Array)
#		 measurements (Array)
#        sig_rep=0 (A scalar that chooses the sigma representation that will be used. By default, its 0)
#	 	 *args will be the Sigma Rep Tunning parameter of sigma representation(May be the weight of sig_rep. Scalar)
#Outputs: x_posterior (Array MxN)
#		  p_posterior (Array NxN)
def ad_ukf(x_previous, pxx_previous, state_funct, param_state_funct, measur_funct, params_measur_funct, Q, R,
		   measurements, sig_rep=0, *args):
	#If there's no sig_rep_tun_param, then it should raise an error
	#*args receives the sigma representation turning parameter
	if (len(args) == 0 and sig_rep != 0):
		raise ValueError("No Sigma Representation Tuning Parameter detected!")

	#This will test the x_previous and measurements array dimensions
	#if they're an array with only one dimension like (1,) instead of (1,N), 
	#then, they will be reshaped to an array Nx1
	dim_state = np.shape(x_previous)
	dim_measur = np.shape(measurements)
	if(len(dim_state) < 2):
		x_previous = np.reshape(x_previous, (dim_state[0], 1))
	if(len(dim_measur) < 2):
		measurements = np.reshape(measurements, (dim_measur[0], 1))


	measur_dimension = measurements.shape[0] # amount of measurements

	#State UT
	x_predicted, p_predicted = ut(x_previous, pxx_previous, state_funct, param_state_funct, args[0], 2, sig_rep)
	p_predicted += Q

	#Measurement UT
	y_predicted, p_measurement, p_cross_covariance = ut(x_predicted, p_predicted, measur_funct, 
														params_measur_funct, args[0], 3, sig_rep, measur_dimension)

	p_measurement += R

	#Kalman Filter correction
	x_posterior, p_posterior = kf_correction(x_predicted, p_predicted, y_predicted, p_measurement,
											 p_cross_covariance, measurements)


	return x_posterior, p_posterior


