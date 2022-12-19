import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv('data/PowerByPlanes.txt', sep=' ')

r=36
Ntheta=200
theta=data.iloc[:,0]
phi=data.iloc[:,1]
P=data.iloc[:,2]
# theta = theta.head(Ntheta)
# phi = phi.head(Ntheta)

x=P*np.sin(theta)*np.cos(phi)
y=P*np.sin(theta)*np.sin(phi)
z=P*np.cos(theta)

fig = plt.figure(figsize = (5,5))
ax = plt.axes(projection="3d")
ax.scatter3D(x,y,z,s=1)
plt.show()
