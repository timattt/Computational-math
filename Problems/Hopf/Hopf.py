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
import matplotlib.cm as cm
from matplotlib.colors import Normalize

def Lagrange(F, X, x):
    n = len(F)
    
    L = 0
    
    for k in range(n):
        P = 1
        for j in range(n):
            if k != j:
                P = P * (x - X[j])/(X[k]-X[j])
        
        L += F[k] * P
        
    return L

def solve(startU, borderU, h, tau, N, M):
    u = [[0 for _ in range(N)] for _ in range(M)]
    xs = np.linspace(0, N*h, N)
    ts = np.linspace(0, M*h, M)
    # H = H[t][x][y]
    
    # start conditions
    
    for i in range(N):
        x = i * h
        u[0][i] = startU(x)
    for n in range(M):
        t = n * tau
        u[n][0] = borderU(t)

    for n in range(M-1):
        for i in range(1, N):
            U = u[n][i]
            for _ in range(3):
                if U > 1:
                    x0 = h - tau / U
                    #S = Lagrange(s[n], xs, x0)
                    #U = Lagrange(u[n], xs, x0)
                    U = Lagrange([u[n][i-1], u[n][i]], [0, h], x0)
                else:
                    t0 = tau - h / U
                    U = Lagrange([u[n][i-1], u[n+1][i-1]], [0, tau], t0)
            
            u[n+1][i] = U
           
    return np.array([np.array([u[n][i] for i in range(N)]) for n in range(M)])

def startU(x):
    if x > 2:
        return 3
    return 5

def borderU(t):
    return 5

h = 0.3
tau = 0.3
N = 100
M = 500
us = solve(startU, borderU, h, tau, N, M)
xs = np.arange(0, N * h, h)
ts = np.arange(0, M * tau, tau)

fig = plt.figure()
ax = fig.add_subplot()

def data_gen(framenumber):
    ax.clear()
    ax.set_ylim(3, 7)
    plots = [ax.plot(xs, us[framenumber % M])]
    return [plots]

plots = [ax.plot(xs, us[0])]
pam_ani = animation.FuncAnimation(fig, data_gen, interval=1, blit=False)
#pam_ani.save("hopf.gif", fps = 10)
plt.show()
