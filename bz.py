import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import odeint
from scipy.integrate import ode
import scipy
from scipy.signal import argrelextrema
eps = [0.11, 0.000017, 0.0016]
gamma = 1.2
alpha = 0.1
beta = 0.000017
mu = 0.00024
k = 0.00037
phi_0 = 0.00016

n = 100               # no. of oscillator
phi_p = 0.147*np.ones(n)
mean = 0.7855
sigma = 0.1086
q = np.random.normal(mean, sigma, n)
# x = [0.98, 0.5]
# y = [0.2, 0.4]
# x = np.array(x)
# y = np.array(y)
# y = q*x

y = np.random.uniform(0,1,n)
x = np.random.uniform(0,1,n)
# y = np.ones(n)*0.6
# x = np.ones(n)
# x = y*q
y0 = np.hstack((x,y))
print(y0)
# def U_ss(X,Z):

# 	Uss = (1/(4*gamma*eps[1]))*( scipy.sqrt(16*gamma*X*eps[1] + Z**2 - 2*Z + 1) + Z - 1 )
# 	print(scipy.sqrt(16*gamma*X*eps[1] + Z**2 - 2*Z + 1))
# 	return Uss 

# Uss = U_ss(1.04, 0.96)
t0 = 500
t_p = 700
def phi_func(t, b, n):
	if ((t > t0) and (t<t0+10)) or ((t>t0+t_p) and (t<t0+t_p+10))  :
		phi = phi_p
		return phi
	else:
		return phi_0 + k*(np.sum(b)/n - b)

def model(t, y, n):
	
	a = y[0:n]
	b = y[n:2*n]
	# print(a[10], b[10])
	phi = phi_func(t, b, n)
	# Uss = U_ss(a,b)
	# print(Uss)
	Uss = (1/(4*gamma*eps[1]))*( np.sqrt(16*gamma*a*eps[1] + b**2 - 2*b + 1) + b - 1 )
	dxdt =(1/eps[0])*(phi - a**2 - a + eps[1]*gamma*Uss**2 + Uss*(1-b) + ((mu-a)/(mu+a))*((q*alpha*b)/(eps[2]+1-b) + beta))
	dzdt =(2*phi + Uss*(1-b) - alpha*b/(eps[2] + 1 - b))

	derv = np.hstack((dxdt,dzdt))
	# derv = np.hstack((dervx,dervz))
	return derv

def theta_func(t, arr):
	peak_pos = t[argrelextrema(arr, np.greater)]
	theta = []
	point = 1
	print(len(peak_pos))
	tn = peak_pos[1]
	tn_1 = peak_pos[0]
	for i in t:
		if (i>peak_pos[1] and i<peak_pos[-1]):
			if(i>peak_pos[point + 1]):
				point = point + 1 
				tn = peak_pos[point]
				tn_1 = peak_pos[point - 1]
				# print(point, i)
			theta.append((i-tn)/(tn - tn_1))
		else:
			theta.append(0)	
	return theta

# print(U_ss(x,z))
t = np.linspace(0, 3000, 10000)
delta_t = t[1] - t[0]
print(delta_t)
arr = []
r = ode(model).set_integrator(name = "lsoda", method = "BDF", nsteps = 100000)
r.set_initial_value(y0, 0).set_f_params(n)
for i in range(len(t)):
	print(i)
	r.integrate(r.t + delta_t)
	arr.append(r.y)
arr = np.array(arr)

thetaArr = []
for k in range(n, 2*n):
	print(k)
	thetaArr.append(theta_func(t, arr[:,k]))
	print(k)
thetaArr = np.array(thetaArr)
R_theta = (1/n)*abs(np.sum(np.exp(2*np.pi*1j*thetaArr), axis = 0))
plt.plot(t, R_theta)
# plt.plot(t, arr[:,n:])
plt.show()
np.savetxt(btp.txt, delimiter=',')