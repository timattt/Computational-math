import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
import math
from Const import *

def task(N, Xmin, Xmax, U, Umin=None, Umax=None, maxDiscrete = None):
    # init values
    h = (Xmax - Xmin) / N
    xs = np.array([Xmin + i *h for i in range(N)])
    U = np.array([U(x) for x in xs])
    if Umin == None:
        Umin = np.min(U)
    if Umax == None:
        Umax = np.max(U)
    solid = -(dirac**2) / (2*m)
    beta = solid / (h**2)
    alpha = np.array([U[i] - 2*solid / (h**2) for i in range(N)])
    A = np.zeros((N, N))
    
    # init matrix
    for i in range(N):
        A[i, i] = alpha[i]
        
    for i in range(N-1):
        A[i, i+1] = A[i+1, i] = beta
    
    # solve spectral problem
    w, v = np.linalg.eigh(A)
    
    print("Umin = {}, Umax = {}".format(Umin, Umax))
    
    # normalize
    psis = np.array([v[:, i] / np.sqrt(np.dot(v[:, i], v[:, i])*h) for i in range(N)])
    psiMax = np.max(np.abs(psis))
    totalDiscrete = np.sum([1 if w[i] < Umax and w[i] > Umin else 0 for i in range(N)])
    
    shown = 0
    
    xs = xs
    for i in range(N):
            psi_ = psis[i]
            E = w[i]
            if E < Umax and E > Umin and (maxDiscrete == None or shown < maxDiscrete):
                psi_ *= (Umax-Umin)/(np.max(psi_)-np.min(psi_))/e/(totalDiscrete if maxDiscrete == None else maxDiscrete)
                psi_ += E/e
                plt.plot(xs/ab, psi_/1000000, label="E={}".format(E/e))
                shown += 1
    
    plt.plot(xs/ab, U/e/1000000)            
    plt.xlabel("радиус бора")
    plt.ylabel("МЭВ")
    print("Total discrete levels: {}".format(totalDiscrete))
    plt.legend()
    plt.show()