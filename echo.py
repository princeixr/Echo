import numpy as np 
import matplotlib.pyplot as plt 
import scipy
from scipy.integrate import quad

X = 0.602
n = 10000     #no. of oscillator
mean = 0.02352941
sigma = 0.001257695
bz_phase = np.random.uniform(0,1, n)
omega    = np.random.normal(mean, sigma, n)
num = np.arange(0,n)

tau = 1000
time = np.arange(7000)
itt = len(time)
lamda_t = np.zeros(itt)
Rt_abs = np.zeros(itt)


def complex_integral(func, a, b, **kwargs):
    def real_func(x):
        return scipy.real(func(x))	
    def imag_func(x):
        return scipy.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1], imag_integral[1])

def const_b(q_value):
	if q_value != 0:
		b_q = (1 - np.exp(-2*np.pi*1j*q_value*X))/(2*np.pi*1j*q_value)
	return b_q

def Coefficient_Cq(q_value, X):
	if q_value == 0:
		return (1-X)*(X - const_b(1))
	if q_value == 1:
		return (X**2)-(1-X)*const_b(1)
	else :
		return X*(const_b(q_value-1) - const_b(q_value)) - (abs(const_b(q_value-1)))**2 + const_b(q_value)*const_b(1-q_value)

C_q = np.zeros(itt//tau+1)
for q in range(itt//tau+1):
	if q==0 :
		C_q[0] = 0
	else:
		C_q[q] = Coefficient_Cq(q-1, X)
print(C_q.shape)
# plt.hist(np.random.normal(0, sigma/8, n))
# plt.show()
G_t_q = np.zeros((itt//tau, itt))

O_p = np.zeros(itt) #Order parameter
R_t = np.zeros(itt)
for t in range(itt):

	q = t//tau	
	omega = omega + np.random.normal(0, sigma, n)*0.006
	bz_phase = bz_phase + omega
	bz_phase[bz_phase > 1] = 0
	bz_phase[bz_phase < 0] = 1
	if ((t>tau) and (t<tau+10)) or ((t>2*tau) and (t<(2*tau + 10))) :
		bz_phase[bz_phase > 1-X] = 0

	else:
		bz_phase = bz_phase
	O_p[t] = (abs(np.sum(np.exp(2*np.pi*1j*bz_phase))))/n

	# if q==0 :
	# 	C_q[q] = 0
	# else:
	# 	C_q[q] = Coefficient_Cq(q-1, X)
	G_t = np.exp(-0.5*(2*np.pi*sigma*(t-q*tau))**2 + 2*np.pi*1j*mean*(t-q*tau))
	R_t[t] = abs(C_q[q])*abs(G_t)

	# real_part = lamda[0].real 
	# imag_part = lamda[0].imag
	# lamda_abs = abs(real_part + 1j*imag_part)
	# # print(lamda_abs)
	# G_t = np.exp(-0.5*(2*np.pi*sigma*(time-q*tau))**2 + 2*np.pi*1j*mean*(time-q*tau))
	# lamda_t[t] = lamda_abs	 
	# if q >= 0 :
	# 	G_t = np.exp(-0.5*(2*np.pi*sigma*(time-q*tau))**2 + 2*np.pi*1j*mean*(time-q*tau))
	# 	C_t = Coefficient_Cq(q,X)
	# 	G_t_q[q:] = G_t
	# 	# print(G_t_q[q:])
	# 	# plt.plot(time, G_t_q[q:].real)	
	# else:
	# 	C_t = 0

plt.plot(omega, bz_phase, 'ro')
plt.show()

plt.plot(time, O_p)
plt.show()



# G_t = [np.exp(-0.5*(2*np.pi*sigma*(time-q*tau))**2 + 2*np.pi*1j*mean*(time-q*tau)) for q in range(itt//tau)]
# G_t = np.array(G_t)
# C_q = np.array(C_q).reshape(itt//tau,1)
# plt.plot(C_q)
# plt.show()

# R_t = C_q*G_t

# R_t_final = np.sum(R_t, axis=0)
# R_t_final = abs(R_t_final)
plt.plot(time, R_t)
plt.show()
# print(R_t_final.shape)










