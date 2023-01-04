import matplotlib.pyplot as plt
import numpy as np
import math

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

def DividedDifference(X, F, k):
    val = 0
    for i in range(k+1):
        vv = F[i]
        for j in range(k+1):
            if j != i:
                vv /= (X[i]-X[j])
        val += vv
    return val

def Newton(X, F):
    n = len(X)
    
    P = np.poly1d([0])

    for i in range(n):
        D = np.poly1d([1])
        for j in range(0, i):
            D = np.polymul(D, np.poly1d([1, -X[j]]))

        P = np.polyadd(P, D * DividedDifference(X, F, i))

    return P.coeffs[::-1]

def Newton_add(X, F, P):
    Q = np.poly1d(P[::-1])
    n = len(X)

    D = np.poly1d([1])
    for j in range(0, n-1):
        D = np.polymul(D, np.poly1d([1, -X[j]]))

    Q = np.polyadd(Q, D * DividedDifference(X, F, n-1))

    return Q.coeffs[::-1]

def Newton_approximate(Func, a, b, n):
    DX = (b - a)/ n
    
    xs = np.arange(a - 0.1, b + 0.1, 0.01)
    fs = []
    fsl = []
    
    X = np.arange(a, b, DX)
    X = X + DX / 2
    
    F = np.array([Func(x) for x in X])
    
    coefs = Newton(X, F)
    
    for x in xs:
        fs.append(Func(x))
        p = 1
        V = 0
        for c in coefs:
            V += c * p
            p *= x
        fsl.append(V)

    fig, ax = plt.subplots()
    
    ax.plot(xs, fsl, label = "interpolated: " + printPol(coefs))
    ax.plot(xs, fs, '--', label = "original")
    ax.scatter(X, F, color = "red", label = "points")

    # add point
    print("Add point X...")
    g = 6.3#float(input())
    X = np.append(X, g)
    F = np.append(F, Func(X[-1]))
    
    coefs = Newton_add(X, F, coefs)

    fsl = []

    for x in xs:
        fs.append(Func(x))
        p = 1
        V = 0
        for c in coefs:
            V += c * p
            p *= x
        fsl.append(V)
        
    ax.plot(xs, fsl, color = "cyan", label = "interpolated +1 : " + printPol(coefs))
    ax.scatter([X[-1]], [F[-1]], color = "green", label = "new point")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    ax.grid()
    ax.legend()
    
    plt.show()
    
def F(X):
    return math.exp(X)*math.sin(X)
    
print("Lagrange approximation function: f(x) = exp(x)sin(x)")
print("Input interval [a, b] and points count N...")
a = 0#float(input())
b = 10#float(input())
n = 5#int(input())

Newton_approximate(F, a, b, n)