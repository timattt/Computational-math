import numpy as np
import matplotlib.pyplot as plt
from Const import *

def task(N, xmin, xmax, U, energy, ns=None, mlt = 1, mayRender = True, mayRound = True, mayDrawEnergyLevels = True):
    print("[ Quasi classics schrodinger solver ]")
    
    if ns == None:
        ns = range(0,5)
    
    print("[ Levels requested: {}".format(len(ns)))
    
    h = (xmax - xmin) / N
    
    xs = np.arange(xmin, xmax, h)
    
    plt.plot(xs/ab, [U(x)/e for x in xs], color = "orange")
    
    du = (np.max( [U(x) for x in xs])-np.min( [U(x) for x in xs])) / len(ns)

    for n in ns:
        En = 0
        x1 = 0
        x2 = 0
            
        En, x1, x2 = energy(n)
    
        dx = (x2-x1)/10000
        
        while En - U(x1) > 0:
            x1 -= dx
        while En - U(x1) <= 0:
            x1 += dx
        while En - U(x2) > 0:
            x2 += dx
        while En - U(x2) <= 0:
            x2 -= dx
        
        print("[ E{}={}, x1={}, x2={}".format(n, round(En/e, 0) if mayRound else En/e, x1/ab, x2/ab))
        
        if x2 < x1:
            print("error!")
            quit(0)
        
        def p(arg):
            if 2*m*(En-U(arg)) < 0:
                print("x={}, U={}".format(arg, U(arg)))
                #quit(0)
            return np.sqrt(2*m*(En-U(arg)))
        
        def I(X):
            res = 0
            x = X
            dx = (x2-X) / N
            while x < x2:
                res += (p(x)*dx)
                x += dx
            return res
        
        def p_comp(arg):
            return np.sqrt(2*m*abs(En-U(arg)))
        
        def I_comp1(X):
            res = 0
            x = X
            dx = (x1-X) / N
            while x < x1:
                res += (p_comp(x)*dx)
                x += dx
            return res
        
        def I_comp2(X):
            res = 0
            x = x2
            dx = (X-x2) / N
            while x < X:
                res += (p_comp(x)*dx)
                x += dx
            return res
                
        def psi(x):
            if x < x1:
                return np.exp(-I_comp1(x)/dirac)/(2*np.sqrt(p_comp(x)))
            if x2 < x:
                return np.exp(-I_comp2(x)/dirac)/(2*np.sqrt(p_comp(x)))
            return np.sin(I(x)/dirac + np.pi/4)/np.sqrt(p(x))
        
        psis = np.array([psi(x) for x in xs])
        
        norm = np.dot(psis, psis) * h
        psis *= 1/np.sqrt(norm)
        psis = psis * psis*du /np.max(psis * psis) * mlt
        psis += En
        
        plt.scatter([x1/ab, x2/ab], [U(x1)/e, U(x2)/e], color="black")
        plt.plot(xs/ab, psis/e, label = "E{}={} ЭВ".format(n, round(En/e, 0) if mayRound else En/e), color = "blue")
        if mayDrawEnergyLevels:
            plt.plot([xs[0]/ab, xs[-1]/ab], [En/e, En/e], linestyle ='--', color="gray")

    plt.xlabel("радиусы бора")
    plt.ylabel("Электронвольты")
    if mayRender:
        plt.legend()
        plt.show()

