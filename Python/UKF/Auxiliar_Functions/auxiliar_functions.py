import numpy as np
import math as mth

#Auxiliar functions implementations: state_or_measure_predicted_cov_matrix_predicted, p_crosscovariance_with_weights, kf_correction, feval
#Author: Victor Araujo Vieira.
#e-mail: icevct@gmail.com
#University of Brasilia - Brazil.


#state_or_measure_predicted_cov_matrix_predicted function definition
#Inputs: sigma_points (an array NxM)
#		 weights (a column array Nx1)
#Outputs: sample_predicted (an array MxN)
#		  p_predicted (an array NxN)
def state_or_measure_predicted_cov_matrix_predicted(sigma_points, weights):
    
    N = weights.shape[0] # amount of sigma points and weights
    #State's or measures prediction
    sample_predicted = np.dot(sigma_points, weights)

    #Covariance matrix prediction
    sig_pts_size = sigma_points.shape[0] #amount of sigma points rows
    p_predicted = np.zeros((sig_pts_size, sig_pts_size))

    for i in range(0, N):
    	p_predicted = p_predicted + weights[[i], :]*np.dot(sigma_points[:, [i]] - sample_predicted, 
    									     np.transpose(sigma_points[:, [i]] - sample_predicted)) 

    return sample_predicted, p_predicted

#p_crosscovariance_with_weights function definition
#Inputs: sigma_points_before (an array XxM)
#		 sigma_points_now (an array YxN)
#		 weights (an row array Nx1)
#Outputs: p_crosscovariance (Array XxY)
def p_crosscovariance_with_weights(sigma_points_before, sigma_points_now, weights):

	#dimensions attributions
	M = sigma_points_now.shape[1] #amount of sigma points_now
	X = sigma_points_before.shape[0] #amount of states
	Y = sigma_points_now.shape[0] # measures dimension

	#state's and measurement prediction
	state_predicted = np.dot(sigma_points_before, weights)
	measurement_predicted = np.dot(sigma_points_now, weights)

	p_crosscovariance = np.zeros((X, Y))

	for i in range(0, M):
		p_crosscovariance = p_crosscovariance + weights[[i], 0]*np.dot(sigma_points_before[:, [i]] - state_predicted,
													  np.transpose(sigma_points_now[:, [i]] - measurement_predicted))

	return p_crosscovariance

#kf_correction function definition
#As there's so many variables, the arrays dimensions were not specified
#Inputs: x_predicted (Array)
#		 p_predicted (Array)
#        y_predicted (Array)
#        p_measurement (Array)
#        p_crosscovariance (Array)
#        measurements (Array)
#Outputs: x_posterior (Array)
#		  p_posterior (Array)
def kf_correction(x_predicted, p_predicted, y_predicted, p_measurement, p_crosscovariance, measurements):

	#Kalman filter correction equation
	Kk = p_crosscovariance/p_measurement
	x_posterior = x_predicted + np.dot(Kk, (measurements - y_predicted))
	p_posterior = p_predicted - np.dot(Kk, np.dot(p_measurement, np.transpose(Kk)))

	return x_posterior, p_posterior

#feval function definition
#N: amount of states
#Inputs: f (A function)
#		 evaluation_values (Values to evaluate the function. An array Nx1)
#		 params (An array 1xN with the parameters)
#Outputs: f_evaluated (An array Nx1 with the function evaluated)
def feval(f, *args):
	if len(args) == 1:
		f_evaluated = f(args[0])
	else:
		f_evaluated = f(args[0], args[1]) 

	return f_evaluated

