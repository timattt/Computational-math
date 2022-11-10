import matplotlib.pyplot as plt
import numpy as np
from Const import *
import Schrodinger

def task(N, Xmin, Xmax, U, V, maxDiscrete = None, mlt = 1, mayRender = True, mayRound = True, mayDrawEnergyLevels = True):
    E, psis, totalDiscrete = Schrodinger.task(N, Xmin, Xmax, U, maxDiscrete, mlt, False, mayRound, mayDrawEnergyLevels=False)
    
    print("[ Stationary perturbation theory solver ]")
    print("[ perturbation percentage={} %".format(round(V(Xmax) / U(Xmax)*100)))
    
    h = (Xmax - Xmin) / N
    
    xs = np.array([Xmin + i *h for i in range(N)])
    
    resultE = [0 for _ in range(N)]
    resultPsis = [0 for _ in range(N)]
    
    Umin = np.min([U(x) for x in xs])
    Umax = np.max([U(x) for x in xs])
    
    du = (Umax - Umin) / (totalDiscrete if maxDiscrete == None else maxDiscrete)

    for n in range(N):
        psi = psis[n]
        En = E[n]
        
        def matrixElem(k, n):
            return np.sum([psis[n][i]*V(x)*psis[k][i]*h for i, x in enumerate(xs)])
        
        # dEn = <n|V|n>
        dEn = matrixElem(n, n)#np.sum([psi[i]*V(x)*psi[i] for i, x in enumerate(xs)])
        
        # |PSIn> = sum (<n|V|k>) / (En - Ek)
        dpsi = np.sum([ (matrixElem(n, k) / (E[n] - E[k]) if k != n else 0) for k in range(N)])
        
        resEn = En + dEn
        resPsi = psi + dpsi
        
        norm = np.sqrt(np.dot(resPsi, resPsi)*h)
        resPsi /= norm
        resPsi = resPsi * resPsi * du / (np.max(resPsi*resPsi)) * mlt
        resPsi += resEn
        
        if resEn < Umax and resEn > Umin and (maxDiscrete == None or n < maxDiscrete):
            print("[ dEn={} ЭВ".format(dEn/e))
            
            plt.plot(xs/ab, resPsi/e, label="E{}={} ЭВ".format(n, round(resEn/e, 0) if mayRound else resEn/e), color="green")
            
            if mayDrawEnergyLevels:
                plt.plot([xs[0]/ab, xs[-1]/ab], [resEn/e, resEn/e], linestyle ='--', color="gray")
        
        resPsi -= resEn
        
        norm = np.sqrt(np.dot(resPsi, resPsi)*h)
        resPsi /= norm
        
        resultE[n] = resEn
        resultPsis[n] = resPsi
 
    plt.plot(xs/ab, np.array([U(x) for x in xs])/e, color = "orange")            
    plt.xlabel("радиусы бора")
    plt.ylabel("Электронвольты")
    if mayRender:
        plt.legend()
        plt.show()
 
    return resultE, resultPsis