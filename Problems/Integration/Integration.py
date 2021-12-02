import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math
import random

def integrateLeftRect(f, a, b, n, draw=True):
    I = 0
    
    if draw:
        fig, ax = plt.subplots()
    
    xs = np.arange(a, b, (b - a) / 1000)
    fs = np.array([f(x) for x in xs])
    
    if draw:
        ax.plot(xs, fs, label = "function", color = "red")
    
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    
    for i in range(n):
        xi = a + i * (b - a) / n
        xip1 = a + (i + 1) * (b - a) / n
        
        dx = xip1 - xi
        I += f(xi) * dx

        rect = patches.Rectangle((xi, 0), dx, f(xi), edgecolor='green', facecolor='none', hatch='/')
        
        if draw:
            ax.add_patch(rect)
        
    if draw:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.legend(["function", "integral: " + str(I)])
        #ax.grid()
        #ax.legend()

        plt.show()
        
    return I

def integrateCenterRect(f, a, b, n, draw=True):
    I = 0
    
    if draw:
        fig, ax = plt.subplots()
    
    xs = np.arange(a, b, (b - a) / 1000)
    fs = np.array([f(x) for x in xs])
    
    if draw:
        ax.plot(xs, fs, label = "function", color = "red")
    
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    
    for i in range(n):
        xi = a + i * (b - a) / n
        xip1 = a + (i + 1) * (b - a) / n
        
        Zetta = (xi+xip1)/2
        
        dx = xip1 - xi
        I += f(Zetta) * dx

        rect = patches.Rectangle((xi, 0), dx, f(Zetta), edgecolor='green', facecolor='none', hatch='/')
        
        if draw:
            ax.add_patch(rect)
        
    if draw:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.legend(["function", "integral: " + str(I)])
        #ax.grid()
        #ax.legend()

        plt.show()
        
    return I

def integrateTrapezoid(f, a, b, n, draw=True):
    I = 0
    
    if draw:
        fig, ax = plt.subplots()
    
    xs = np.arange(a, b, (b - a) / 1000)
    fs = np.array([f(x) for x in xs])
    
    if draw:
        ax.plot(xs, fs, label = "function", color = "red")
    
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    
    for i in range(n):
        xi = a + i * (b - a) / n
        xip1 = a + (i + 1) * (b - a) / n
        
        fxi = f(xi)
        fxip1 = f(xip1)
        
        dx = xip1 - xi
        I += 0.5 * (fxi + fxip1) * dx

        if draw:
            ax.add_patch(patches.Polygon(xy=list(zip( [xi, xip1, xip1, xi]  , [0, 0, fxip1, fxi] )), fill=False, color = "green", hatch='/'))
        
    if draw:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.legend(["function", "integral: " + str(I)])
        #ax.grid()
        #ax.legend()
        
        plt.show()
        
    return I

def integrateSimpson(f, a, b, n, draw=True):
    I = 0
    
    if draw:
        fig, ax = plt.subplots()
    
    xs = np.arange(a, b, (b - a) / 1000)
    fs = np.array([f(x) for x in xs])
    
    if draw:
        ax.plot(xs, fs, label = "function", color = "red")
    
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    
    for i in range(n):
        xi = a + i * (b - a) / n
        xip1 = a + (i + 1) * (b - a) / n
        xiph = (xi+xip1)/2
        
        fxi = f(xi)
        fxip1 = f(xip1)
        fxiph = f(xiph)
        
        dx = xip1 - xi
        I += (1/6) * (fxi + 4 * fxiph + fxip1) * dx

        xs = np.arange(xi, xip1, 0.01)
        ys = np.array([fxi*(x-xiph)*(x-xip1)/((xi-xiph)*(xi-xip1)) + fxiph * (x-xi)*(x-xip1)/((xiph-xi)*(xiph-xip1)) + fxip1*(x-xi)*(x-xiph)/((xip1-xi)*(xip1-xiph)) for x in xs])
        
        if draw:
            ax.plot(xs, ys, color = "green")
            ax.stem([xi, xip1], [fxi, fxip1])
            ax.scatter([xiph], [fxiph], color="blue")
    
    if draw:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.legend(["function", "integral: " + str(I)])
        #ax.grid()
        #ax.legend()
        plt.show()
        
    return I

def integrateMontecarlo(f, a, b, n, draw=True):
    I = 0
    
    if draw:
        fig, ax = plt.subplots()
    
    xs = np.arange(a, b, (b - a) / 1000)
    fs = np.array([f(x) for x in xs])
    
    y0 = min(fs)
    y1 = max(fs)
    
    if draw:
        ax.plot(xs, fs, label = "function", color = "green")
    
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    
    k = 0
    
    xsr = []
    ysr = []
    
    xsb = []
    ysb = []
    
    for i in range(n):
        x = random.uniform(a, b)
        y = random.uniform(y0, y1)
        
        v = f(x)
       
        if (v >= y and v >= 0 and y >= 0) or (v < 0 and y >= v and y < 0):
            xsr.append(x)
            ysr.append(y)
            k += 1
        else:
            xsb.append(x)
            ysb.append(y)    
              
    I = (b - a) * (y1 - y0) * k / n + y0*(b-a)
                
    if draw:
        ax.scatter(xsb, ysb, color="blue")
        ax.scatter(xsr, ysr, color="red", label="integral: " + str(I))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        
        #ax.grid()
        ax.legend()
        
        plt.show()
        
    return I

def testRungeRule(func, trueInt, a, b, N):
    fig, ax = plt.subplots()
    
    errs1 = []
    errs2 = []
    errs3 = []
    errs4 = []
    
    rerrs1 = []
    rerrs2 = []
    rerrs3 = []
    rerrs4 = []
    
    start = 4
    
    for n in range(start, N):
        I1 = integrateLeftRect(func, a, b, 2*n, False)
        I2 = integrateCenterRect(func, a, b, 2*n, False)
        I3 = integrateTrapezoid(func, a, b, 2*n, False)
        I4 = integrateSimpson(func, a, b, 2*n, False)
        
        I1_ = integrateLeftRect(func, a, b, n, False)
        I2_ = integrateCenterRect(func, a, b, n, False)
        I3_ = integrateTrapezoid(func, a, b, n, False)
        I4_ = integrateSimpson(func, a, b, n, False)
        
        err1 = math.fabs(trueInt - I1)
        err2 = math.fabs(trueInt - I2)
        err3 = math.fabs(trueInt - I3)
        err4 = math.fabs(trueInt - I4)
        
        rerr1 = (1/3)*math.fabs(I1-I1_)
        rerr2 = (1/3)*math.fabs(I2-I2_)
        rerr3 = (1/3)*math.fabs(I3-I3_)
        rerr4 = (1/15)*math.fabs(I4-I4_)
        
        errs1.append(err1)
        errs2.append(err2)
        errs3.append(err3)
        errs4.append(err4)
        
        rerrs1.append(rerr1)
        rerrs2.append(rerr2)
        rerrs3.append(rerr3)
        rerrs4.append(rerr4)
        
    ns = [2*n for n in range(start, N)]
    
    ax.plot(ns, errs1, label = "left rect", color = "green")
    ax.plot(ns, errs2, label = "Center rect", color = "blue")
    ax.plot(ns, errs3, label = "Trapezoid", color = "red")
    ax.plot(ns, errs4, label = "Simpson", color = "orange")
    
    ax.plot(ns, rerrs1, label = "Runge feft rect", linestyle='dashed', color = "green")
    ax.plot(ns, rerrs2, label = "Runge center rect", linestyle='dashed', color = "blue")
    ax.plot(ns, rerrs3, label = "Runge trapezoid", linestyle='dashed', color = "red")
    ax.plot(ns, rerrs4, label = "Runge simpson", linestyle='dashed', color = "orange")
    
    ax.set_xlabel("n")
    ax.set_ylabel("error")
    
    ax.legend()
    plt.show()

Z = 2

def func(x):
    return math.sin(x)+Z

a = 0
b = 10
N = 4

integrateLeftRect(func, a, b, N)
integrateCenterRect(func, a, b, N)
integrateTrapezoid(func, a, b, N)
integrateSimpson(func, a, b, N)
integrateMontecarlo(func, a, b, N*N*N)

testRungeRule(func, (math.cos(a) - math.cos(b)+Z*(b-a)), a, b, 100)
