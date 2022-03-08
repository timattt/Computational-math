import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.patches as patches
import numpy as np
import math
import random
import scipy.optimize

def sweep(A, b):
    n = len(b)
    y = [0 for row in range(n)]
    alpha = [0 for row in range(n)]
    beta = [0 for row in range(n)]
    
    if A[0][0] == 0:
        return 0
    
    y[0] = A[0][0]
    alpha[0] = -A[0][1] / y[0]
    beta[0] = b[0] / y[0]
    
    # forward
    for i in range(1, n):
        y[i] = A[i][i] + alpha[i-1] * A[i][i-1]
        
        if y[i] == 0:
            return 0
        
        if i < n - 1:
            alpha[i] = -A[i][i+1] / y[i]
        beta[i] = (b[i] - A[i][i-1]*beta[i-1])/y[i]
        
    # backward and result
    b[n-1] = beta[n-1]
    
    for i in range(n-2, -1, -1):
        b[i] = alpha[i] * b[i+1] + beta[i]
    return 1

def gridMethodLinear(a, b, YaTrue, YbTrue, h, p, q, f):
    N = int((b - a) / h)
    
    xs = np.arange(a, b, h)
    xs = np.append(xs, 1)
    A = [1/h**2 + p(xs[i]) / (2*h)  for i in range(N)]
    B = [-2/h**2 + q(xs[i]) for i in range(N)]
    C = [1/h**2 - p(xs[i]) / (2*h) for i in range(N)]
    D = [f(xs[i]) for i in range(N)]
    
    MAT = []
    VEC = []
    
    MAT.append([0 for i in range(N+1)])
    MAT[0][0] = 1
    VEC.append(YaTrue)
    
    for i in range(1, N):
        ADD = [0 for i in range(N+1)]
        ADD[i - 1] = C[i]
        ADD[i] = B[i]
        ADD[i + 1] = A[i]
        MAT.append(ADD)
        VEC.append(D[i])
        
    MAT.append([0 for i in range(N+1)])
    MAT[-1][-1] = 1
    VEC.append(YbTrue)
    
    sweep(MAT, VEC)
    
    ys = np.array(VEC)
    
    return xs, ys


def FVA(k, V, f, h):
    N = int(1 / h)
    
    xs = np.arange(0, 1, h)
    xs = np.append(xs, 1)

    def V_(i, half):
        if half == 0.5:
            return V((xs[i + 1] + xs[i])/2) 
        if half == 0:
            return V(xs[i])
        if half == -0.5:
            return V((xs[i-1] + xs[i])/2)
    
    A = [(k/h**2 + V_(i, -0.5)/(2*h)) for i in range(N)]
    B = [-( 2*k/h**2 - V_(i, 0.5)/(2*h) + V_(i, -0.5)/(2*h)) for i in range(N)]
    C = [-(V_(i, 0.5) / (2 * h) - k/h**2)  for i in range(N)]
    D = [f(xs[i]) for i in range(N)]
    
    MAT = []
    VEC = []
    
    MAT.append([0 for i in range(N+1)])
    MAT[0][0] = 1
    VEC.append(0)
    
    for i in range(1, N):
        ADD = [0 for i in range(N+1)]
        ADD[i - 1] = C[i]
        ADD[i] = B[i]
        ADD[i + 1] = A[i]
        MAT.append(ADD)
        VEC.append(D[i])
        
    MAT.append([0 for i in range(N+1)])
    MAT[-1][-1] = 1
    VEC.append(0)
    
    sweep(MAT, VEC)
    
    ys = np.array(VEC)

    return xs, ys, [V_(i, 0) * h / k for i in range(0, N+1)]

def monotonizedScheme(k, V, f, h):
    N = int(1 / h)
    
    xs = np.arange(0, 1, h)
    xs = np.append(xs, 1)

    def V_(i, half):
        if half == 0.5:
            return V((xs[i + 1] + xs[i])/2) 
        if half == 0:
            return V(xs[i])
        if half == -0.5:
            return V((xs[i-1] + xs[i])/2)
        
    def k_(i):
        ro = max(0, V(xs[i]) * h / (2 * k) - 1 + 0.01)
        return k * (1 + ro)
    
    A = [(k_(i)/h**2 + V_(i, -0.5)/(2*h)) for i in range(N)]
    B = [-( 2*k_(i)/h**2 - V_(i, 0.5)/(2*h) + V_(i, -0.5)/(2*h)) for i in range(N)]
    C = [-(V_(i, 0.5) / (2 * h) - k_(i)/h**2)  for i in range(N)]
    D = [f(xs[i]) for i in range(N)]
    
    Pe = np.array([V_(i, 0) * h / k_(i) for i in range(0, N+1)])
    
    MAT = []
    VEC = []
    
    MAT.append([0 for i in range(N+1)])
    MAT[0][0] = 1
    VEC.append(0)
    
    for i in range(1, N):
        ADD = [0 for i in range(N+1)]
        ADD[i - 1] = C[i]
        ADD[i] = B[i]
        ADD[i + 1] = A[i]
        MAT.append(ADD)
        VEC.append(D[i])
        
    MAT.append([0 for i in range(N+1)])
    MAT[-1][-1] = 1
    VEC.append(0)
    
    sweep(MAT, VEC)
    
    ys = np.array(VEC)

    return xs, ys, Pe

def DMS(k, V, f, h):
    N = int(1 / h)
    
    xs = np.arange(0, 1, h)
    xs = np.append(xs, 1)

    def V_(i, half):
        if half == 0.5:
            return V((xs[i + 1] + xs[i])/2) 
        if half == 0:
            return V(xs[i])
        if half == -0.5:
            return V((xs[i-1] + xs[i])/2)
    
    def VP(i, half):
        return 0.5 * (V_(i, half) + math.fabs(V_(i, half)))
    
    def VM(i, half):
        return 0.5 * (V_(i, half) - math.fabs(V_(i, half)))
    
    A = [(k/h**2 + VM(i, 0.5)/h) for i in range(N)]
    B = [(2*k/h**2 + VP(i, 0.5)/h - VM(i, -0.5)/h) for i in range(N)]
    C = [(-VP(i, -0.5) / h - k/h**2)  for i in range(N)]
    D = [f(xs[i]) for i in range(N)]
    
    MAT = []
    VEC = []
    
    MAT.append([0 for i in range(N+1)])
    MAT[0][0] = 1
    VEC.append(0)
    
    for i in range(1, N):
        ADD = [0 for i in range(N+1)]
        ADD[i - 1] = C[i]
        ADD[i] = B[i]
        ADD[i + 1] = A[i]
        MAT.append(ADD)
        VEC.append(D[i])
        
    MAT.append([0 for i in range(N+1)])
    MAT[-1][-1] = 1
    VEC.append(0)
    
    sweep(MAT, VEC)
    
    ys = np.array(VEC)

    return xs, ys, [V_(i, 0) * h / k for i in range(0, N+1)]


def V(x):
    return math.exp(10*x)

def f(x):
    return math.exp(-x)

def makeGraph(func, h):
    xs, ys, Pe = func(1, V, f, h)


    fig, axs = plt.subplots(2)

    axs[1].plot(xs, Pe, label = "Peclet")
    i0 = 0
    for i in range(len(xs)):
        if (Pe[i] > 2):
            i0 = i
            break
    
    axs[1].scatter(xs[i0], Pe[i0])

    axs[0].scatter(xs[i0], ys[i0])    

    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")

    axs[1].set_xlabel("x")
    axs[1].set_ylabel("Pe")

    axs[0].plot(xs, ys, label = "solution")

    xs, ys = gridMethodLinear(0, 1, 0, 0, h, V, lambda x: 0, f)

    #axs[0].plot(xs, ys, '--', label = "true")

    axs[0].legend()
    axs[1].legend()
    plt.show()

def makeGraph2(func, h):
    xs, ys, Pe = func(1, V, f, h)


    fig, axs = plt.subplots(2)

    axs[1].plot(xs, Pe, label = "Peclet")
    i0 = 0
    for i in range(len(xs)):
        if (Pe[i] > 2):
            i0 = i
            break
    
    #axs[1].scatter(xs[i0], Pe[i0])

    #axs[0].scatter(xs[i0], ys[i0])    

    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")

    axs[1].set_xlabel("x")
    axs[1].set_ylabel("Pe")

    axs[0].plot(xs, ys, label = "solution")

    xs, ys, tmp = FVA(1, V, f, h)

    axs[0].plot(xs, ys, '--', label = "true")

    axs[0].legend()
    axs[1].legend()
    plt.show()
    
def makeGraph3(func, h):
    xs, ys, Pe = func(1, V, f, h)


    fig, ax = plt.subplots(1)


    ax.set_xlabel("x")
    ax.set_ylabel("y")


    ax.plot(xs, ys, label = "solution")

    xs, ys = gridMethodLinear(0, 1, 0, 0, h, V, lambda x: 0, f)

    #axs[0].plot(xs, ys, '--', label = "true")

    ax.legend()

    plt.show()

#makeGraph(FVA, 0.001)
#makeGraph2(monotonizedScheme, 0.1)
makeGraph3(DMS, 0.05)
