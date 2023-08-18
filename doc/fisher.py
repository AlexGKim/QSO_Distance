import numpy
import scipy.integrate as integrate

flat = False

Omega_M=0.3
Omega_L = (1-Omega_M)

def E(z):
	return numpy.sqrt(Omega_M*(1+z)**3 + Omega_L)

def invE(z):
	return 1/E(z)

def partial_OM_int(z):
	# if flat:
	# 	return ((1+z)**3-1)/E(z)**3
	# else:
	return ((1+z)**3)/E(z)**3

def partial_OL_int(z):
	return 1/E(z)**3



def fisher(zs, sigs):

	coeff = -2.5/numpy.log(10)

	F = numpy.zeros((3,3))

	for z,sig in zip(zs,sigs):
		factor = integrate.quad(invE, 0, z)[0]
		_OM = integrate.quad(partial_OM_int, 0, z)[0]
		_OL = integrate.quad(partial_OL_int, 0, z)[0]		
		F[0,0]  = F[0,0]+ 1./sig**2
		F[0,1]  = F[0,1]+ coeff/factor/sig**2*(_OM + 1./3*factor**3)
		F[0,2]  = F[0,2]+ coeff/factor/sig**2*(_OL + 1./3*factor**3)
		F[1,1]  = F[1,1]+ coeff**2 /factor**2 /sig**2 * (_OM + 1./3*factor**3)**2
		F[1,2]  = F[1,2]+ coeff**2 /factor**2 /sig**2 * (_OM + 1./3*factor**3) *(_OL + 1./3*factor**3)
		F[2,2] = F[2,2]+ coeff**2 /factor**2 /sig**2 * (_OL + 1./3*factor**3)**2

	F[1,0]=F[0,1]
	F[2,0]=F[0,2]
	F[2,1]=F[1,2]
	return F


if __name__ == "__main__":

	# Quasar
	zs = numpy.linspace(.5,4,10)
	sigs = numpy.zeros(len(zs))+0.1

	# A supernova-like example
	# zs = numpy.linspace(0.1,.5,10)
	# sigs = numpy.zeros(len(zs))+0.1
	# zs = numpy.append(zs,0.001)
	# sigs = numpy.append(sigs,0.01)


	f = fisher(zs,sigs)

	f_inv=numpy.linalg.inv(f)
	# print(f)
	# print(numpy.sqrt(1/f.diagonal()))

	print(f_inv)
	print(numpy.sqrt(f_inv.diagonal()))
	print("Uncertainty in Omega_M {}".format(numpy.sqrt(f_inv[1,1])))
	print("Uncertainty in Omega_L {}".format(numpy.sqrt(f_inv[2,2])))