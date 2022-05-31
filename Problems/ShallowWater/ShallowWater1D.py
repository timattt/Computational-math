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
    
    if x+0.01 < X[0] or x > X[-1]:
        print(str(X[0]) + " " + str(x) + " " + str(X[-1]))
        g = x  / 0
        print("err")
    
    L = 0
    
    for k in range(n):
        P = 1
        for j in range(n):
            if k != j:
                P = P * (x - X[j])/(X[k]-X[j])
        
        L += F[k] * P
        
    return L

def solve(startU, startS, borderU, borderS, h, tau, N, M):
    s = [[0 for _ in range(N)] for _ in range(M)]
    u = [[0 for _ in range(N)] for _ in range(M)]
    Rs = [[[0 for _ in range(N)] for _ in range(M)] for _ in range(2)]
    xs = np.linspace(0, N*h, N)
    ts = np.linspace(0, M*h, M)
    # H = H[t][x][y]
    
    ro = 1
    def funcC(s):
        return math.sqrt(math.sqrt(math.fabs(s)) / ro)
    
    # start conditions
    
    for i in range(N):
        x = i * h
        S = s[0][i] = startS(x)
        U = u[0][i] = startU(x)
        C = funcC(S)
        l1 = [-C, S]
        l2 = [C, S]
        Rs[0][0][i] = l1[0]*S + l1[1]*U
        Rs[1][0][i] = l2[0]*S + l2[1]*U
        
    for n in range(M):
        t = n * tau
        S = s[n][0] = s[n][-1] = borderS(t)
        U = u[n][0] = u[n][-1] = borderU(t)
        C = funcC(S)
        l1 = [-C, S]
        l2 = [C, S]
        Rs[0][n][0] = Rs[0][n][-1] = l1[0]*S + l1[1]*U
        Rs[1][n][0] = Rs[1][n][-1] = l2[0]*S + l2[1]*U
    
    for n in range(M-1):
        for i in range(1, N-1):
            R1 = Rs[0][n][i]
            R2 = Rs[1][n][i]
            S = 0
            C = 0
            U = 0
            lambda1 = 0
            
            for k in range(3):
                S = math.pow(math.pow((R1-R2)/2, 4), 0.2) * math.sqrt(ro)
                C = funcC(S)

                # if math.fabs(S) < 0.0001:
                #     R1 = Rs[0][n+1][i-1]
                #     R2 = Rs[1][n+1][i-1]
                #     #U = u[n+1][i-1]
                #     continue
                
                U = (R1 + R2) / (2*S)
                
                lambda1 = U + C
                lambda2 = U - C

                # R1
                if lambda1 <= 0:
                    R1 = R1
                    if lambda1 < -tau / (h * (N-i)):
                        R1 = R1
                        #R1 = Lagrange(Rs[0][n], xs, i * h - tau / lambda1)
                    else:
                        t0 = lambda1 * (N * h - i * h) + tau * (n+1)
                        #R1 = Lagrange([Rs[0][n][-1], Rs[0][n+1][-1]], [n*tau , tau*(n+1)], t0)
                else:
                    if lambda1 > 1:
                        x0 = h - tau / lambda1
                        R1 = Lagrange([Rs[0][n][i-1], Rs[0][n][i]], [0, h], x0)
                    else:
                        t0 = tau - h * lambda1
                        R1 = Lagrange([Rs[0][n][i-1], Rs[0][n+1][i-1]], [0, tau], t0)
                        
                # R2
                if lambda2 <= 0:
                    R2 = R2
                    if lambda2 < -tau / (h * (N-i)):
                        R2 = R2
                        #R2 = Lagrange(Rs[1][n], xs, i * h - tau / lambda2)
                    else:
                        t0 = lambda2 * (N * h - i * h) + tau * (n+1)
                        #R2 = Lagrange([Rs[1][n][-1], Rs[1][n+1][-1]], [n*tau , tau*(n+1)], t0)
                else:
                    if lambda2 > 1:
                        x0 = h - tau / lambda2
                        R2 = Lagrange([Rs[1][n][i-1], Rs[1][n][i]], [0, h], x0)
                    else:
                        t0 = tau - h * lambda2
                        R2 = Lagrange([Rs[1][n][i-1], Rs[1][n+1][i-1]], [0, tau], t0)
                    
                
            Rs[0][n+1][i] = R1
            Rs[1][n+1][i] = R2
            s[n+1][i] = S
            u[n+1][i] = U

    return np.array([np.array([s[n][i] for i in range(N)]) for n in range(M)]), np.array([np.array([u[n][i] for i in range(N)]) for n in range(M)])

def startU(x):
    return 2*math.exp(-x)

def startS(x):
    return 4

def borderS(t):
    return 4

def borderU(t):
    return 2

h = 0.1
tau = 0.1
N = 100
M = 500
ss, us = solve(startU, startS, borderU, borderS, h, tau, N, M)
xs = np.linspace(0, h * N, N)
ts = np.linspace(0, M * tau, M)

fig = plt.figure()
ax = fig.add_subplot()

def data_gen(framenumber):
    ax.clear()
    #ax.set_ylim(-6, 6)
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
#pam_ani.save("pipe.gif")
plt.show()
