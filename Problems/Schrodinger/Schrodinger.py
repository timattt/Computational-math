import matplotlib.pyplot as plt
import numpy as np
from Const import *

def task(N, Xmin, Xmax, U, maxDiscrete = None, mlt = 1, mayRender = True, mayRound = True, mayDrawEnergyLevels = True):
    print("[ Schrodinger spectral solver ]")
    
    # init values
    h = (Xmax - Xmin) / N
    xs = np.array([Xmin + i *h for i in range(N)])
    U = np.array([U(x) for x in xs])
    Umin = np.min(U)
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

    # normalize
    psis = np.array([v[:, i] / np.sqrt(np.dot(v[:, i], v[:, i])*h) for i in range(N)])
    totalDiscrete = np.sum([1 if w[i] < Umax and w[i] > Umin else 0 for i in range(N)])

    print("[ Total discrete levels: {}".format(totalDiscrete))
    
    n = 0
    du = (Umax - Umin) / (totalDiscrete if maxDiscrete == None else maxDiscrete)

    for i in range(N):
            psi_ = psis[i]
            E = w[i]
            if E < Umax and E > Umin and (maxDiscrete == None or n < maxDiscrete):
                print("[ E{}={} ЭВ".format(n, round(E/e, 0) if mayRound else E/e))
                psi_ = psi_ * psi_ * du / (np.max(psi_*psi_)) * mlt
                psi_ += E
                plt.plot(xs/ab, psi_/e, label="E{}={} ЭВ".format(n, round(E/e, 0) if mayRound else E/e), color="red")
                if mayDrawEnergyLevels:
                    plt.plot([xs[0]/ab, xs[-1]/ab], [E/e, E/e], linestyle ='--', color="gray")
                psi_ -= E
                psi_ /= np.sqrt(np.dot(psi_, psi_)*h)
                    
                n += 1
    
    plt.plot(xs/ab, U/e, color = "orange")            
    plt.xlabel("радиусы бора")
    plt.ylabel("Электронвольты")
    if mayRender:
        plt.legend()
        plt.show()
        
    return w, psis, totalDiscrete