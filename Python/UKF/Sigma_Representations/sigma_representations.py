import numpy as np
import math as mth
import scipy.linalg

#Sigma representation implementations: even_homi_sy_sig_rep, ho_mi_sy_sig_rep, rho_mi_sig_rep and sig_rep_julier_1995
#Author: Victor Araujo Vieira.
#e-mail: icevct@gmail.com
#University of Brasilia - Brazil.

#All the implementations are acording to the paper 2015_Menegaz et al._A Systematization of the Unscented Kalman Filter Theory
#Especificaly, the table I

#even_homi_sy_sig_rep function definition
#This is the default sigma representation
#Inputs: state (a column array Mx1)
#		 cov_matrix (Array NxN)
#Outputs: sigma_points_before (Array Mx2*M)
#		  weights (a column array Nx1)
def even_homi_sy_sig_rep(state, cov_matrix):

	N = state.shape[0] # amount of states
	sigma_points_before = np.zeros((N, 2*N))

	delta_sigma_point = np.linalg.cholesky(N*cov_matrix).T #generates the delta_sigma_point
	weights = np.tile(1./(2.*N), (2*N, 1)) #generates the new weights

	#sigma points generation
	for i in range(0, N):
		sigma_points_before[:, [i]] = state + delta_sigma_point[:, [i]]
		sigma_points_before[:, [i + N]] = state - delta_sigma_point[:, [i]]	

	return sigma_points_before, weights

#ho_mi_sy_sig_rep function definition
#(Table I, 2)
#Inputs: state (a column array Mx1)
#		 cov_matrix (Array NxN)
#	 	 weight (scalar)
#Outputs: sigma_points_before (Array Mx(2*M + 1))
#		  weights (a column array (2*M + 1)x1)
def ho_mi_sy_sig_rep(state, cov_matrix, weight):

	if weight > 1:
		raise ValueError('Weight must have a value lesser than 1!')

	n = state.shape[0] #amount of states
	N = 2*n + 1 # amount of sigma points

	weights = np.zeros((N, 1))
	sigma_points_before = np.zeros((n, N))

	#weights generation
	for i in range(0, N):
		weights[[i], 0] = (1 - weight)/(2*n)

	weights[0, 0] = weight #first weight must be the weight passed 
						   #as argument to the function

	delta_sigma_point = np.linalg.cholesky((n/(1 - weight))*cov_matrix).T #generates the delta_sigma_point
	#sigma points generation
	sigma_points_before[:, [0]] = state # first sigma point must be the state
	for i in range(0, n):
		sigma_points_before[:, [i + 1]] = state + delta_sigma_point[:, [i]]
		sigma_points_before[:, [i + n + 1]] = state - delta_sigma_point[:, [i]]

	return sigma_points_before, weights

#rho_mi_sig_rep function definition
#(Table I, 4)
#Inputs: state (a column array Mx1)
#		 cov_matrix (Array NxN)
#	 	 weight (scalar)
#Outputs: sigma_points_before (Array Mx(M + 1))
#		  weights (a column array (M + 1)x1)
def rho_mi_sig_rep(state, cov_matrix, weight):

	if weight < 0 or weight > 1:
		raise ValueError('Weight must have a value higher than 0 or lesser than 1!')

	n = state.shape[0] #amount of states
	N = n + 1 #number of sigma points

	weights = np.zeros((N, 1))
	rho_matrix = np.ones((n, n))
	one_vector = np.ones((n, 1))
	identity_matrix = np.eye(n)
	sigma_points_before = np.zeros((n, N))

	rho = mth.sqrt((1 - weight)/n) # rho definition
	rho_matrix = (rho**2)*rho_matrix # rho_matrix multiplied by rho
	c = np.linalg.cholesky(identity_matrix - rho_matrix).T.T
	aux_matr = np.dot(np.linalg.inv(c), np.dot(weight*rho_matrix, 
		              np.linalg.inv(np.transpose(c))))

	#weights generation
	for i in range(0, n):
		weights[[i], 0] = aux_matr[[i], [i]]

	w = np.diag(weights[0:n, 0])
	weights[N-1, 0] = weight # the last weight must be the weight
	                         #passed as argument to the function

	#auxiliar matrices generation
	aux_matr1 = np.dot(np.linalg.cholesky(cov_matrix).T, 
							  np.dot(c, np.linalg.inv(np.linalg.cholesky(w).T)))
	aux_matr2 = np.dot((-rho)*np.linalg.cholesky(cov_matrix).T, 
		               one_vector/mth.sqrt(weight))
	final_matrix = np.concatenate((aux_matr1, aux_matr2), axis = 1) # pode dar erro, TODO

	#sigma points generation
	for i in range(0, N):
		sigma_points_before[:, [i]] = final_matrix[:, [i]] + state
	
	return sigma_points_before, weights

#sig_rep_julier_1995 function definition
#(Table I, 1)
#Inputs: state (a column array Mx1)
#		 cov_matrix (Array NxN)
#	 	 kappa (scalar)
#Outputs: sigma_points_before (Array Mx(2*M + 1))
#		  weights (a column array (2*M + 1)x1)
def sig_rep_julier_1995(state, cov_matrix, kappa):

	n = state.shape[0] #amount of states
	N = 2*n + 1 #amount of sigma points
	weights = np.zeros((N, 1))
	sigma_points_before = np.zeros((n, N))

	if kappa <= -n:
		raise ValueError('Weight must have a value higher than -n!')

	#weight generation
	new_weight = kappa/(n + kappa)
	for i in range(0, N):
		weights[[i], 0] = 1/(2*(n + kappa))

	weights [0, 0] = new_weight #first weight must be the new_weight calculated
	
	#sigma points generation
	delta_sigma_point = np.linalg.cholesky((n + kappa)*cov_matrix).T #generates the delta_sigma_point

	sigma_points_before[:, [0]] = state # the first sigma point must be the state
	for i in range(0, n):
		sigma_points_before[:, [i + 1]] = state + delta_sigma_point[:, [i]]
		sigma_points_before[:, [i + n + 1]] = state - delta_sigma_point[:, [i]]

	return sigma_points_before, weights