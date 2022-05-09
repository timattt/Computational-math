import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import math
import random
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from mayavi import mlab
from numpy import linalg as LA
from matplotlib import animation


def solve(startH, startU, startV, h, tau, N, M):
    g = 1
    H = [[[0 for _ in range(N+M)] for _ in range(N+M)] for _ in range(M)]
    u = [[[0 for _ in range(N+M)] for _ in range(N+M)] for _ in range(M)]
    v = [[[0 for _ in range(N+M)] for _ in range(N+M)] for _ in range(M)]
    # H = H[t][x][y]
    
    # start conditions
    
    for i in range(N+M):
        for k in range(N+M):
            x = i * h
            y = k * h
            H[0][i][k] = startH(x, y)
            u[0][i][k] = startU(x, y)
            v[0][i][k] = startV(x, y)
    
    HH = 1
    f = 1
    g = 1
    b = 1
    C = tau/h
    for n in range(M-1):
        for i in range(N+M-1-n):
            for k in range(N+M-1-n): 
                #H[n+1][i][k] = H[n][i][k] - C * (H[n][i+1][k]*u[n][i+1][k] - H[n][i][k]*u[n][i][k]) - C * (H[n][i][k+1]*v[n][i][k+1]-H[n][i][k]*v[n][i][k])
                #u[n+1][i][k] = H[n][i][k]*u[n][i][k] - C * (H[n][i+1][k]*u[n][i+1][k]**2 + 0.5*g*H[n][i+1][k]**2 - H[n][i][k]*u[n][i][k]**2 - 0.5*g*H[n][i][k]**2) - C*(H[n][i][k+1]*u[n][i][k+1]*v[n][i][k+1]-H[n][i][k]*u[n][i][k]*v[n][i][k])
                #u[n+1][i][k] /= H[n+1][i][k]
                #v[n+1][i][k] = H[n][i][k]*v[n][i][k] - C*(H[n][i+1][k]*u[n][i+1][k]*v[n][i+1][k] - H[n][i][k]*u[n][i][k]*v[n][i][k]) - C * (H[n][i][k+1]*v[n][i][k+1]**2 + 0.5*g*H[n][i][k+1]**2 - H[n][i][k]*v[n][i][k]**2 - 0.5*g*H[n][i][k]**2)
                #v[n+1][i][k] /= H[n+1][i][k]
                H[n+1][i][k] = H[n][i][k] - HH*C * (u[n][i+1][k] - u[n][i][k] + v[n][i][k+1]-v[n][i][k])
                u[n+1][i][k] = tau*f*v[n][i][k] + u[n][i][k] - g * C * (H[n][i+1][k]-H[n][i][k])-b*tau*u[n][i][k]
                v[n+1][i][k] = -tau*f*u[n][i][k] - g * C * (H[n][i][k+1]-H[n][i][k])-b*tau*v[n][i][k]
    return np.array([np.array([np.array([H[n][i][k] for i in range(N)]) for k in range(N)]) for n in range(M)])

def startH(x, y):
    return 0.1*math.exp(-2*(x-0.5)**2-2*(y-0.5)**2)

def startV(x, y):
    return 0

h = 0.1
tau = 0.01
N = 10
M = 200
Hs = solve(startH, startV, startV, h, tau, N, M)
xs = np.arange(0, N * h, h)
ys = np.arange(0, N * h, h)
ts = np.arange(0, M * tau, tau)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(xs, ys)
   
def data_gen(framenumber):
    ax.clear()
    plot = ax.plot_surface(X, Y, Hs[framenumber%M])
    return [plot,]
    
plot = ax.plot_surface(X, Y, Hs[0])
pam_ani = animation.FuncAnimation(fig, data_gen, interval=30, blit=False)
plt.show()
