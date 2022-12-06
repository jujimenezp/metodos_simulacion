import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

epsilon0 = 1.
mu0 = 2.
T=25
C = 1/np.sqrt(2)
alpha = 0.25
J0=1

Z0 = np.sqrt(mu0/epsilon0)
omega = 2*np.pi/T
k = omega/C
s = 1/np.sqrt(2*alpha)
p = J0/omega*(s**3)*(2*np.pi)**1.5

A = Z0*k*k*p/(4*np.pi)

By = lambda x,t: A*abs(np.cos(k*x-omega*t) + np.sin(k*x-omega*t)/(k*x))/x

x = np.arange(58,100,0.1)
y = By(x,70)

fig, ax = plt.subplots()
plt.rcParams['figure.figsize']=(10,5)
ax.plot(x,y)
ax.grid()

fig.savefig("perfil.png")
plt.show()
