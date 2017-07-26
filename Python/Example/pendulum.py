import numpy as np
import math as mth
import scipy.linalg
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from Ad_UKF import *

#PENDULUM EXAMPLE

#This code was adapted from MATLAB to python by me.
#The original one was coded by Simo Sarkka

# Simulate pendulum data for the examples in the book
#
# Simo Sarkka (2013), Bayesian Filtering and Smoothing,
# Cambridge University Press. 
#
#.

# Simulate simple pendulum. Note that the system easily
# diverges, but it should not matter.

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



def pendulum():
	DT = 0.01
	g = 9.81
	Q = 0.01*np.array([[(DT**3)/3, (DT**2)/2], [(DT**2)/2, DT]])
	R = 0.1
	m0 = np.array([[1.6], [0]])
	P0 = 0.1*np.eye(2)

	#defining the state function parameters
	param_F = np.zeros((2))
	param_F[0] = DT
	param_F[1] = g

	steps = 500

	QL = np.linalg.cholesky(Q)

	y = np.zeros((1, steps))
	X = np.zeros((2, steps))
	x = np.array([[1.5], [0]])
	x2 = np.zeros((2, 1))
	w = np.zeros((2, 2))
	T = np.zeros((1, steps))
	t = 0
	x2 = x

	#data generation
	for k in range(0, steps):	
		x[0, 0] = x2[0, 0] + DT*x2[1, 0]
		x[1, 0] = x2[1, 0] - g*mth.sin(x2[0, 0])*DT
		w = np.dot(QL, np.random.randn(2, 1))
		x += w
		x2 = x
		X[:, [k]] = x
		t += DT
		y[0, [k]] = mth.sin(x[0, 0]) + mth.sqrt(R)*np.random.randn(1, 1)
		T[0, [k]] = t


	#plotting the data generated 
	plt.plot(T[0, :], y[0, :], 'r.', label = 'Measurements')
	plt.plot(T[0, :], X[0, :], 'b-', label = 'True angle')
	plt.legend()
	plt.title('Measurements and True angle')	
	#plt.plot(x[0, :])
	plt.show()

	x_prev = m0
	pxx_prev = P0
	state = np.zeros((2, 1))
	MM = np.zeros((x_prev.shape[0], steps))

	#ukf computation
	for k in range(0, steps):
		state, p = ad_ukf(x_prev, pxx_prev, state_function, param_F, measurement_function, [], Q, R,
					  y[0, [k]], 3, 0.1)
		x_prev = state
		pxx_prev = p
		MM[:, [k]] = state


	plt.plot(T[0, :], y[0, :], 'r.', label = 'Measurements')
	plt.plot(T[0, :], X[0, :], 'b-', label = 'True angle')
	plt.plot(T[0, :], MM[0, :], 'g-', label = 'UKF')
	plt.legend()
	plt.title('After UKF')
	#plt.text(1, -1.9,'Red - Measurements, Blue - Real State, Green - Estimated State')
	plt.show()


pendulum()