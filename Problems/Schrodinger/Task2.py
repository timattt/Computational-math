from Const import ab, dirac, m, e
import numpy as np
import Schrodinger
import QuasiClassics
import matplotlib.pyplot as plt
import Task2Tmp

N = 700

def packedTask(Xmin, Xmax, a, U0, ns):
    Task2Tmp.first = Task2Tmp.second = Task2Tmp.third = 0
    
    print("[ alpha0={}".format(np.sqrt(2*m*U0)/dirac))
    print("[ U0={}".format(U0/e))
    
    def U(x):
        return U0 * (a/x - x/a)**2
        
    def energy(n):
        import Task2Tmp
        
        alpha0 = np.sqrt(2*m*U0)/dirac

        p1 = alpha0/10
        p2 = alpha0*10
        
        En = -1
        x1 = -1
        x2 = -1
        
        if n < p1:
            Task2Tmp.first += 1
            En = 2 * np.sqrt(U0)*dirac*(n+0.5)/ (a*np.sqrt(2*m))
            x1 = a*(1-np.sqrt(En/(4*U0)))
            x2 = a*(1+np.sqrt(En/(4*U0)))
        if n > p1 and n < p2:
            Task2Tmp.second += 1
            En = 8 * dirac**2 *(n+0.5)**2 / (m * a**2)
            x1 = a*(np.sqrt(5)/2 - 0.5)
            x2 = a*(np.sqrt(5)/2 + 0.5)
        if n > p2:
            Task2Tmp.third += 1
            En = 4 * np.sqrt(U0)*dirac*(n+0.5)/ (a*np.sqrt(2*m))
            x1 = np.sqrt(U0/En)*a
            x2 = np.sqrt(En/U0)*a
            
        return En, x1, x2
    
    
    #Schrodinger.task(N, Xmin, Xmax, U, maxDiscrete=4, mlt=0.2, mayRender=False, mayDrawEnergyLevels=False)
    QuasiClassics.task(N, Xmin, Xmax, U, energy, ns, mlt=1, mayRender=False, mayRound=False)
    
    print("[ n << alpha0 -> {}, n ~ alpha0 -> {}, n >> alpha0 -> {}\n".format( Task2Tmp.first, Task2Tmp.second, Task2Tmp.third))
    
    plt.legend()
    plt.show()
    
    
packedTask(0.45, 0.55, 0.5, 0.00000000000001*e, range(1,5)) 
packedTask(0.06, 3.5, 0.5, 0.000000000000000001*e, range(1,5))
packedTask(0.018, 15, 0.5, 0.0000000000000000000001*e, range(1, 5))
packedTask(0.04, 6, 0.5, 0.000000000000002*e, [1, 2, 100, 200, 3000, 4000]) 