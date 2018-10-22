import numpy as np
from matplotlib import pyplot as plt

def phi (theta,beta):
	FHS=2.545
	FHS=27.20614
	FCHS=1
	ratio=FHS/FCHS
	K=ratio**(-1/3.7)
	return theta - np.arccos((1-K*(1+beta*np.cos(theta*2*np.pi/360)))/beta)*360/(2*np.pi);

#def deltatrig(theta):
#	return np.arcsin(np.sin(theta*np.pi/180)/4)*180/(np.pi);

def deltatrig(theta):
	return np.arcsin(np.sin(theta*np.pi/180)/2)*180/(np.pi);

def phi_d(theta):
	return theta-deltatrig(theta);
minb = 0.4185
maxb = 0.4205
brange = np.linspace(minb,maxb,num=int(np.abs((maxb-minb)/0.00002)))
theta1 = np.linspace(-1,5.5,num=1000)

def phi_d_max(theta):
	return phi_d(theta) + theta;


for i in range(len(brange)):
	if i==0:
		plt.plot(theta1,phi(theta1,brange[i]),c='Blue',label=r'$\varphi$',linewidth = 0.1)
	else:
		plt.plot(theta1,phi(theta1,brange[i]),c='Blue',linewidth = 0.1)

plt.plot(theta1,theta1,c='Orange', label='bisector')
#plt.plot(theta1,deltatrig(theta1),c='Red', label=r'$\delta^{\prime}$')
plt.plot(theta1,phi_d(theta1),c='Green', label=r'$\varphi^{\prime}$')
plt.plot(theta1,phi_d_max(theta1),'--',c="Green", label=r'$\varphi^{\prime}_{max}$')
plt.fill_between(theta1,phi_d(theta1),phi_d_max(theta1),facecolor='Green',alpha=0.2)
plt.ylim(0,6)
plt.xlim(0,5.5)
plt.xlabel(r'$\theta$ [deg]')
plt.ylabel(r'respective angle [deg]')
plt.legend()


plt.savefig('calculusparam_max.pdf')
plt.show()
