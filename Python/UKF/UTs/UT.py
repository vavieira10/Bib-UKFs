import numpy as np
import math as mth
import sys 
sys.path.append('..')
from Auxiliar_Functions.auxiliar_functions import *
from Sigma_Representations.sigma_representations import *

#Unscented transformation implementation: ut
#Author: Victor Araujo Vieira.
#e-mail: icevct@gmail.com
#University of Brasilia - Brazil.

#ut function definition
#As there's so many variables, the arrays dimensions were not specified
#Inputs: x_previous (Array)
#		 pxx_previous (Array)
#        f (A function. The state or measurement state function)
#        param_f (Array or an empty list [])
#        weight (scalar)
#		 nargout (scalar): the variable that controls how many 
#							values the function will return
#        *args (scalar): It can be 1 or 2 argument, which will be the option to which sigma represetantion will be used or 
#						 dim_y a scalar, that indicates the amount of measurements 
#Outputs: x_predicted (Array)
#		  p_predicted (Array)
#		  p_cross_covariance (Array)
def ut(x_previous, pxx_previous, f, param_f, weight, nargout, *args):
	
	if(len(args) == 1):
		n = x_previous.shape[0]; # amount of states
	else:
		n = args[1] #amount of measurements

	#Chooses which sigma representation will be used, based on the option passed as argument
	if(args[0] == 0):
	#even_homi_sy_sig_rep
	    sigma_points_before, weights = even_homi_sy_sig_rep(x_previous, pxx_previous)
	elif(args[0] == 1):
	#sig_rep_julier_1995
	    sigma_points_before, weights = sig_rep_julier_1995(x_previous, pxx_previous, weight)
	elif(args[0] == 2):
	#ho_mi_sy_sig_rep
	    sigma_points_before, weights = ho_mi_sy_sig_rep(x_previous, pxx_previous, weight)
	elif(args[0] == 3):
    #rho_mi_sig_rep
	    sigma_points_before, weights = rho_mi_sig_rep(x_previous, pxx_previous, weight)
	#elif(args[0] == 4):
	else:
	    raise ValueError("There's not that option for a Sigma Representation!")


	N = sigma_points_before.shape[1] #amount of weights and sigma points columns 
	sigma_points_now = np.zeros((n, N)) 
	#sigma points of current sample, evaluating the sigma points to the state or measurement function

	for i in range(0, N):
		if(param_f == []):
			sigma_points_now[:, [i]] = feval(f, sigma_points_before[:, [i]])
		else:
			sigma_points_now[:, [i]] = feval(f, sigma_points_before[:, [i]], param_f)

	#function that calculates the x_predicted and p_predicted
	x_predicted, p_predicted = state_or_measure_predicted_cov_matrix_predicted(sigma_points_now, weights)

	#function that generates the cross covariance matrix
	p_cross_covariance = p_crosscovariance_with_weights(sigma_points_before, sigma_points_now, weights)

	#if nargout = 2, that means that there's only 2 return values
	#if nargout = 3, that means that there's 3 return values
	#the AdUKF function (main) will deal with it
	if nargout == 2:
		return x_predicted, p_predicted
	elif nargout == 3:
		return x_predicted, p_predicted, p_cross_covariance

