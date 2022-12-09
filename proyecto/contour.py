import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('data/antena_70.dat', sep='\t')
x = data.iloc[:,0]
y = data.iloc[:,1]
By = data.iloc[:,2]

fig,ax=plt.subplots(1,1)
contours = ax.tricontour(x,y,By,colors='black',linewidths=0.6)
By_array = np.resize(By,(100,100))
cp=ax.imshow(By_array,extent=[0,100,0,100])
fig.colorbar(cp)
ax.set_xlim([5,95])
ax.set_ylim([5,95])
plt.show()
