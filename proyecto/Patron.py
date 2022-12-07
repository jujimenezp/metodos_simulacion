import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

dataH = pd.read_csv('data/patronH.dat', sep='\t')
dataE = pd.read_csv('data/patronE.dat', sep='\t')

thetaH = dataH.iloc[:,0]
PH = dataH.iloc[:,1]

thetaE = dataE.iloc[:,0]
PE = dataE.iloc[:,1]

theta2 = np.arange(0,2*np.pi,0.01)
r=-6.4*np.ones(theta2.size)

epsilon0 = 1.; mu0 = 2.;
lamb=16; C = 1/np.sqrt(2); T=lamb/C; omega = 2*np.pi/T
alpha = 0.5; J0=0.0001; Z0 = np.sqrt(mu0/epsilon0); s = 1/np.sqrt(2*alpha)
I=J0*(np.pi/alpha)
dP=0.4*np.log(Z0*I*I/(8*np.pi**2)*(np.cos(np.pi/2*np.sin(theta2))**2)/(np.cos(theta2)**2))

fig,axs=plt.subplots(1,2,subplot_kw={'projection': 'polar'})
axs[0].plot(thetaH,PH,'rx')
axs[0].plot(theta2,r,'b')
axs[0].grid(True)
axs[0].set_rticks([0,-1,-2,-3,-4,-5,-6,-7])
axs[0].set_ylim([0,-8])

axs[1].plot(thetaE,PE,'rx')
axs[1].plot(theta2,dP,'b')
axs[1].grid(True)
#axs[1].set_rticks([0,-10,-20,-30,-40])
#axs[1].set_ylim([-50,0])
plt.show()
