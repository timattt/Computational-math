import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def prob(r, te):
    E = 1
    return (np.exp(-r)/np.sqrt(np.pi) - E * r * (2+r)*np.exp(-r)/np.sqrt(4*np.pi) * np.cos(te))**2

x = y = np.arange(0, 3, 0.01)
X, Y = np.meshgrid(x, y)
zs = np.array([[prob(r, te) for te in y] for r in x])
Z = zs.reshape(X.shape)

sur = ax.plot_surface(X, Y, Z,cmap=cm.coolwarm)

ax.set_xlabel('r')
ax.set_ylabel('tetta')
ax.set_zlabel('probability')

plt.show()
    