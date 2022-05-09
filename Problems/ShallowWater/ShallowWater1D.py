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

def solve(startU, startS, borderU, borderS, funcC, h, tau, N, M):
    s = [[0 for _ in range(N)] for _ in range(M)]
    u = [[0 for _ in range(N)] for _ in range(M)]
    xs = np.linspace(0, N*h, N)
    ts = np.linspace(0, M*h, M)
    # H = H[t][x][y]
    
    # start conditions
    
    for i in range(N):
        x = i * h
        s[0][i] = startS(x)
        u[0][i] = startU(x)
    for n in range(M):
        t = n * tau
        s[n][0] = borderS(t)
        u[n][0] = borderU(t)

    for n in range(M-1):
        for i in range(1, N):
            S = s[n][i]
            U = u[n][i]
            C = funcC(S)
            for k in range(3):
                x0 = h - tau * C / S
                if x0 > 0:
                    x0 += i * h - h
                    #S = Lagrange(s[n], xs, x0)
                    #U = Lagrange(u[n], xs, x0)
                    S = Lagrange([s[n][i-1], s[n][i]], [i * h - h, i * h], x0)
                    U = Lagrange([u[n][i-1], u[n][i]], [i * h - h, i * h], x0)
                else:
                    S = Lagrange([s[n][i-1], s[n+1][i-1]], [0, tau], -x0 * S/C)
                    U = Lagrange([u[n][i-1], u[n+1][i-1]], [0, tau], -x0 * S/C)
                    
                C = funcC(S)    
            
            s[n+1][i] = S
            u[n+1][i] = U
           
    return np.array([np.array([s[n][i] for i in range(N)]) for n in range(M)]), np.array([np.array([u[n][i] for i in range(N)]) for n in range(M)])

def startU(x):
    return 0

def startS(x):
    return 2

def funcC(s):
    return math.sqrt(math.sqrt(math.fabs(s)))

def borderS(t):
    return 4+1.5*math.sin(t)

def borderU(t):
    return 4+1.5*math.sin(t)

h = 0.1
tau = 0.5
N = 50
M = 500
ss, us = solve(startU, startS, borderU, borderS, funcC, h, tau, N, M)
xs = np.arange(0, N * h, h)
ts = np.arange(0, M * tau, tau)

fig = plt.figure()
ax = fig.add_subplot()

def data_gen(framenumber):
    ax.clear()
    ax.set_ylim(-6, 6)
    plots = [ax.plot(xs, ss[framenumber % M]), ax.plot(xs, -ss[framenumber % M])]
    mx = max(np.fabs(us[framenumber % M]))
    if mx == 0:
        mx = 1
    for i in range(1, len(xs)):
        x1 = xs[i-1]
        x2 = xs[i]
        s1 = ss[framenumber%M][i-1]
        s2 = ss[framenumber%M][i]
        plots.append(ax.fill([x1, x2, x2, x1], [s1, s2, -s2, -s1], color = (0, 0, np.fabs(us[framenumber % M][i]) / mx)))
    return [plots]

plots = [ax.plot(xs, ss[0]), ax.plot(xs, -ss[0])]
pam_ani = animation.FuncAnimation(fig, data_gen, interval=1, blit=False)
pam_ani.save("pipe.gif")
plt.show()
