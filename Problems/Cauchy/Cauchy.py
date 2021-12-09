import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.patches as patches
import numpy as np
import math
import random
import scipy.optimize

def EulerExplicit(func, x0, y0, n, b):
    x = np.arange(x0, b, (b - x0) / n)
    y = [0 for _ in x]
    y[0] = y0
    
    for i in range(1, len(x)):
        y[i] = y[i - 1] + (x[i] - x[i-1])*func(x[i-1], y[i-1])
        
    return x, y

def EulerImplicit(func, x0, y0, n, b):
    x = np.arange(x0, b, (b - x0) / n)
    y = [0 for _ in x]
    y[0] = y0
    
    for i in range(1, len(x)):
        y[i] = scipy.optimize.fsolve(lambda yi: yi - y[i-1] - (x[i] - x[i-1])*func(x[i-1], yi), y[i-1]) 
        
    return x, y

def EulerRecount(func, x0, y0, n, b):
    x = np.arange(x0, b, (b - x0) / n)
    y = [0 for _ in x]
    y[0] = y0
    
    for i in range(1, len(x)):
        y_ = y[i - 1] + (x[i] - x[i-1])*func(x[i-1], y[i-1])
        y[i] = y[i-1] + (x[i] - x[i-1])*(func(x[i-1], y[i-1]) + func(x[i], y_))*0.5
        
    return x, y

def RungeKutta(f, x0, y0, n, b):
    h = (b - x0) / n
    x = np.arange(x0, b, h)
    y = [0 for _ in x]
    y[0] = y0
    
    for i in range(1, len(x)):
        k1 = f(x[i-1], y[i-1])
        k2 = f(x[i-1] + h/2, y[i-1] + 0.5*h*k1)
        k3 = f(x[i-1] + h/2, y[i-1] + 0.5*h*k2)
        k4 = f(x[i-1] + h, y[i-1] + h * k3)
        
        y[i] = y[i-1] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        
    return x, y

def testErrors(func, trueFunc):
    x0 = 0
    y0 = 1
    b = 1
    
    def err(n, method):
        x, y = method(func, x0, y0, n, b)
        er = 0
        for i in range(len(x)):
            er = max(er, math.fabs(y[i] - trueFunc(x[i])))
        return er
    
    errsEE = []
    errsEI = []
    errsER = []
    errsRK = []
    
    ns = [n for n in range(3, 400)]
    
    for n in ns:
        errsEE.append(err(n, EulerExplicit))
        errsEI.append(err(n, EulerImplicit))
        errsER.append(err(n, EulerRecount))
        errsRK.append(err(n, RungeKutta))
    
    fig, ax = plt.subplots()
    ax.set_xlabel("n")
    ax.set_ylabel("err")

    ax.plot(ns, errsEE, color = "blue", label = "Euler explicit")
    ax.plot(ns, errsEI, color = "green", label = "Euler implicit")
    ax.plot(ns, errsER, color = "red", label = "Euler recount")
    ax.plot(ns, errsRK, color = "yellow", label = "Runge-Kutta")
    
    ax.grid()
    ax.legend()
    ax.set_title("Errors")
        
    plt.show()
    
    ns = [math.log(n) for n in ns]
    errsEE = [math.log(er) for er in errsEE]
    errsEI = [math.log(er) for er in errsEI]
    errsER = [math.log(er) for er in errsER]
    errsRK = [math.log(er) for er in errsRK]
    
    fig, ax = plt.subplots()
    ax.set_xlabel("ln(n)")
    ax.set_ylabel("ln(err)")

    ax.plot(ns, errsEE, color = "blue", label = "Euler explicit ~ O(h)")
    ax.plot(ns, errsEI, color = "green", label = "Euler implicit ~ O(h)")
    ax.plot(ns, errsER, color = "red", label = "Euler recount ~ O(h^2)")
    ax.plot(ns, errsRK, color = "yellow", label = "Runge-Kutta ~ O(h^4)")
    
    ax.grid()
    ax.legend()
    ax.set_title("Log errors")
        
    plt.show()
    

def testSolution(func, equation, trueFunc, method, methodName, x0, y0, n, b, ystep, N):
    
    fig, ax = plt.subplots()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    color = iter(cm.rainbow(np.linspace(0, 1, N)))
    
    for i in range(N):
        x, y = method(func, x0, y0, n, b)
        c = next(color)
        if N < 25 or i % (N / 20) == 0:
            ax.plot(x, y, color = c, label = "f(" + str(x0) + ") = " + str(y0))
        else:
            ax.plot(x, y, color = c)
        
        y0 += ystep
        
    if trueFunc != None:
        xs = [x for x in np.arange(x0, b, 0.001)]
        ys = [trueFunc(x) for x in xs]
        ax.plot(xs, ys, color = "black", label = "true solution for y(0) = 1")
        
    ax.set_title(equation + "\n" + methodName)
    ax.grid()
    ax.legend()
        
    plt.show()

def func0(x, y):
    return x*y

def trueFunc0(x):
    return math.exp(x*x/2)

def func1(x, y):
    return x*math.sin(x)*y

def trueFunc1(x):
    return math.exp(math.sin(x) - x * math.cos(x))

def func2(x, y):
    return math.sin(x*y*y)

def func3(x, y):
    return math.sin(x)*math.exp(-y)

def trueFunc3(x):
    return math.log(2 - math.cos(x))

def func4(x, y):
    return y**math.sin(x)

testSolution(func1, "dy/dx = x*y*sin(x)", trueFunc1, EulerExplicit, "Euler explicit", 0, 1, 1000, 5, 0.01, 200)
testSolution(func2, "dy/dx = sin(x*y*y)", None, EulerImplicit, "Euler implicit", 0, 1, 1000, 1.5, 0.001, 100)
testSolution(func3, "dy/dx = sin(x)*exp(-y)", trueFunc3, EulerRecount, "Euler recount", 0, 0, 1000, 10, 0.01, 300)
testSolution(func4, "dy/dx = y^x", None, RungeKutta, "Runge-Kutta", 0, 0, 1000, 20, 0.005, 200)

testErrors(func0, trueFunc0)