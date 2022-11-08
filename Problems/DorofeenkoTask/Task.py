import math
import numpy as np
import matplotlib.pyplot as plt

dirac = 1#(10**(-34))
m = 1#9.1 * (10**(-31))
ab = 1#0.5 * (10**(-10))
e = 1#1.6 * (10**(-19))

N = 1000

def task(xmin, xmax, a, U0, ns=None):
    h = (xmax - xmin) / N
    
    def U(x):
        return U0 * (a/x - x/a)**2
    
    xs = np.arange(xmin, xmax, h)
    
    alpha0 = np.sqrt(2*m*U0)/dirac
    
    p1 = alpha0/10
    p2 = alpha0*10
    
    print("alpha0={}".format(alpha0))
    
    plt.plot(xs, [U(x) for x in xs], label = "U(x)")
    
    du = (np.max( [U(x) for x in xs])-np.min( [U(x) for x in xs])) / 10
    
    first = 0
    second = 0
    third = 0
    
    if ns == None:
        ns = range(1,5)
    
    for n in ns:
        En = 0
        x1 = 0
        x2 = 0
        if n < p1:
            first += 1
            En = 2 * np.sqrt(U0)*dirac*(n+0.5)/ (a*np.sqrt(2*m))
            x1 = a*(1-np.sqrt(En/(4*U0)))
            x2 = a*(1+np.sqrt(En/(4*U0)))
        if n > p1 and n < p2:
            second += 1
            En = 8 * dirac**2 *(n+0.5)**2 / (m * a**2)
            x1 = a*(np.sqrt(5)/2 - 0.5)
            x2 = a*(np.sqrt(5)/2 + 0.5)
        if n > p2:
            third += 1
            En = 4 * np.sqrt(U0)*dirac*(n+0.5)/ (a*np.sqrt(2*m))
            x1 = np.sqrt(U0/En)*a
            x2 = np.sqrt(En/U0)*a
    
        dx = (x2-x1)/10000
        
        while En - U(x1) > 0:
            x1 -= dx
        while En - U(x1) <= 0:
            x1 += dx
        while En - U(x2) > 0:
            x2 += dx
        while En - U(x2) <= 0:
            x2 -= dx
        
        print("En={}, x1={}, x2={}".format(En, x1, x2))
        
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
            dx = (x2-x1) / N
            while x < x2:
                res += (p(x)*dx)
                x += dx
            return res
        
        def p_comp(arg):
            return np.sqrt(2*m*abs(En-U(arg)))
        
        def I_comp1(X):
            res = 0
            x = X
            dx = (x2-x1) / N
            while x < x1:
                res += (p_comp(x)*dx)
                x += dx
            return res
        
        def I_comp2(X):
            res = 0
            x = x2
            dx = (x2-x1) / N
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
        
        psis = [psi(x)**2 for x in xs]
        
        norm = np.dot(psis, psis) * h
        psis /= np.max(np.abs(psis))#np.sqrt(norm)
        psis *= du
        psis += En
        
        #plt.plot(xs, [En for x in xs])
        plt.scatter([x1, x2], [U(x1), U(x2)], color="black")
        plt.xlabel("радиус бора")
        plt.ylabel("электрон-вольт")
        plt.plot(xs, psis, label = "E={}".format(En))
        
        
    print("n << alpha0 -> {}, n ~ alpha0 -> {}, n >> alpha0 -> {}\n".format( first, second, third))
    plt.legend()
    plt.show()
   
task(4.4*ab, 5.7*ab, 5*ab, 1000*e) 
task(1*ab, 17*ab, 5*ab, 1*e)
task(1.5*ab, 70*ab, 10*ab, 0.001*e)
task(1.3*ab, 17*ab, 5*ab, 100*e, [1, 10, 144])