import matplotlib.pyplot as plt
import numpy as np
import math

def drawGraph(coefs, X, F):
    xs = np.arange(X[0] - 0.1, X[-1] + 0.1, 0.01)
    fs = []

    for x in xs:
        val = 0
        p = 1
        for i in range(len(coefs)):
            val += p*coefs[i]
            p *= x
    
        fs.append(val)
    
    fig, ax = plt.subplots()
    
    ax.plot(xs, fs, label = "interpolated")
    ax.scatter(X, F, color = "red", label = "original")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    ax.grid()
    ax.legend()
    
    plt.show()

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

def Lagrange_coef(X, F):
    n = len(F)
    
    L = []
    
    for k in range(n):
        P = np.array([0 for i in range(2*n)])
        P[0] = 1
        P[1] = 0
        for j in range(n):
            if k != j:
                a = 1/(X[k]-X[j])
                b = -X[j]/(X[k]-X[j])
                
                P = np.polymul(P, [b, a]);
        L = np.polyadd(L, P * F[k])
        
    return L

def printPol(coefs):
    items = []
    for i, x in enumerate((coefs)):
        if not x:
            continue
        items.append('{}x^{}'.format("{0:0.2f}".format(x) if x != 1 or i == 0 else '', i))
    result = ' + '.join(items)
    result = result.replace('x^0', '')
    result = result.replace('^1 ', ' ')
    result = result.replace('+ -', '- ')
    
    return result

def Lagrange_approximate(Func, a, b, n):
    DX = (b - a)/ n
    
    xs = np.arange(a - 0.1, b + 0.1, 0.01)
    fs = []
    fsl = []
    
    X = np.arange(a, b, DX)
    X = X + DX / 2
    F = np.array([Func(x) for x in X])
    
    for x in xs:
        fs.append(Func(x))
        fsl.append(Lagrange(F, X, x))

    fig, ax = plt.subplots()
    
    ax.plot(xs, fsl, label = "interpolated: " + printPol(Lagrange_coef(X, F)))
    ax.plot(xs, fs, '--', label = "original")
    ax.scatter(X, F, color = "red", label = "points")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    ax.grid()
    ax.legend()
    
    plt.show()
    
    printPol(Lagrange_coef(X, F))

def F(X):
    return math.exp(X)*math.sin(X)
    
print("Lagrange approximation function: f(x) = exp(x)sin(x)")
print("Input interval [a, b] and points count N...")
a = float(input())
b = float(input())
n = int(input())

Lagrange_approximate(F, a, b, n)