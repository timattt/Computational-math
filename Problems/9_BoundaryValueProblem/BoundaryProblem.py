import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.patches as patches
import numpy as np
import math
import random
import scipy.optimize

# OLD
#=======================================================
def differentiate_2or_1der(points:list, step:float) -> list:
    if len(points) < 3:
        raise RuntimeError("too short array of points!")
    
    res = []
    n = len(points)
    
    # from left to right
    for i in range(n - 2):
        f0 = points[i]
        f1 = points[i+1]
        f2 = points[i+2]
        der = ( -1.5*f0 + 2*f1 - 0.5*f2 ) / step
        res.append(der);
    
    # right border
    for i in range(n-2, n, 1):
        f0 = points[i]
        f_1 = points[i-1]
        f_2 = points[i-2]
        der = (0.5 * f_2 - 2 * f_1 + 1.5 * f0) / step
        res.append(der)
   
    return res

def RungeKutta(f, x0, y0, h, b):
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

def newton_chords(x0, x1, epsilon, f):
    xn = x0
    xnm1 = x1
    
    for i in range(1000000):
        val_fxn = f(xn)
        val_fxn1 = f(xnm1)
        xnp1 = xn - val_fxn * (xn - xnm1) / (val_fxn - val_fxn1)
        if math.fabs(xnp1 - xn) < epsilon:
            break
        xnm1 = xn
        xn = xnp1
    
    return xn

def bin_search(a, b, epsilon, f):
    while b - a > epsilon:
        c = (a + b) / 2
        if f(b) * f(c) < 0:
            a = c
        else:
            b = c
    f((a + b) / 2)
    return (a + b) / 2

def sweep(A, b):
    n = len(b)
    y = [0 for row in range(n)]
    alpha = [0 for row in range(n)]
    beta = [0 for row in range(n)]
    
    if A[0][0] == 0:
        return 0
    
    y[0] = A[0][0]
    alpha[0] = -A[0][1] / y[0]
    beta[0] = b[0] / y[0]
    
    # forward
    for i in range(1, n):
        y[i] = A[i][i] + alpha[i-1] * A[i][i-1]
        
        if y[i] == 0:
            return 0
        
        if i < n - 1:
            alpha[i] = -A[i][i+1] / y[i]
        beta[i] = (b[i] - A[i][i-1]*beta[i-1])/y[i]
        
    # backward and result
    b[n-1] = beta[n-1]
    
    for i in range(n-2, -1, -1):
        b[i] = alpha[i] * b[i+1] + beta[i]
    return 1


#=======================================================

def shootingMethodGeneral(a, b, YaTrue, YbTrue, h, Y_, Y__):
    YS = []
    XS = []   
    
    def Yb(p):
        xs, ys = RungeKutta(lambda X, Y: np.array([Y_(X, Y[0], Y[1]), Y__(X, Y[0], Y[1])]), a, np.array([YaTrue, p]), h, b)
        XS.append(xs)
        YS.append(ys)
        return (ys[-1][0] - YbTrue)
    
    Ya_ = bin_search(-100, 100, 0.000001, Yb)#newton_chords(-100, 100, h, Yb)
    #print(Ya_)
    return XS, YS, Ya_

def shootingMethodLinear(a, b, YaTrue, YbTrue, h, p, q, f):
    G = 0
    xs, ys_homogeneous = RungeKutta(lambda X, Y: np.array([Y[1], -p(X) * Y[1] + -q(X) * Y[0]]), a, np.array([YaTrue, 2]), h, b) 
    xs, ys_private = RungeKutta(lambda X, Y: np.array([Y[1], -p(X) * Y[1] + -q(X) * Y[0] + f(X)]), a, np.array([YaTrue, 4]), h, b) 

    y_ha = ys_homogeneous[0][0]
    y_hb = ys_homogeneous[-1][0]
    
    y_pa = ys_private[0][0]
    y_pb = ys_private[-1][0]

    #alpha, beta = np.linalg.solve(np.array([[y_ha, y_pa], [y_hb, y_pb]]), np.array([YaTrue, YbTrue]))
    
    alpha = (YbTrue - y_pb) / y_hb
    
    ys = [alpha * ys_homogeneous[i] + ys_private[i] for i in range(len(xs))]
    #print(ys[0][1])
    return [xs, xs, xs], [ys_homogeneous, ys_private, ys], ys[0][1]

def gridMethodLinear(a, b, YaTrue, YbTrue, h, p, q, f):
    N = int((b - a) / h)
    
    xs = np.arange(a, b, h)
    A = [1/h**2 + p(xs[i]) / (2*h)  for i in range(N)]
    B = [-2/h**2 + q(xs[i]) for i in range(N)]
    C = [1/h**2 - p(xs[i]) / (2*h) for i in range(N)]
    D = [f(xs[i]) for i in range(N)]
    
    MAT = []
    VEC = []
    
    MAT.append([0 for i in range(N+1)])
    MAT[0][0] = 1
    VEC.append(YaTrue)
    
    for i in range(1, N):
        ADD = [0 for i in range(N+1)]
        ADD[i - 1] = C[i]
        ADD[i] = B[i]
        ADD[i + 1] = A[i]
        MAT.append(ADD)
        VEC.append(D[i])
        
    MAT.append([0 for i in range(N+1)])
    MAT[-1][-1] = 1
    VEC.append(YbTrue)
    
    sweep(MAT, VEC)
    
    ys = np.array(VEC)
    ys_der = differentiate_2or_1der(ys, h)
    #print(ys_der[0])
    return [xs], [[np.array([ys[i], ys_der[i]]) for i in range(len(xs))]], ys_der[0]

def fixedPointsMethodLinear(a, b, YaTrue, YbTrue, h, p, q, f):
    N = int((b - a) / h)
    
    xs = np.arange(a, b + h, h)
    A = [1/h**2 + p(xs[i]) / (2*h)  for i in range(N+1)]
    B = [-2/h**2 + q(xs[i]) for i in range(N+1)]
    C = [1/h**2 - p(xs[i]) / (2*h) for i in range(N+1)]
    D = [f(xs[i]) for i in range(N+1)]
    
    MAT = []
    VEC = []
    
    MAT.append([0 for i in range(N+2)])
    MAT[0][0] = -1 / (2*h)
    MAT[0][1] = 1 / (2*h)
    VEC.append(YaTrue)
    
    for i in range(0, N):
        ADD = [0 for i in range(N+2)]
        ADD[i - 1 + 1] = C[i]
        ADD[i + 1] = B[i]
        ADD[i + 1 + 1] = A[i]
        MAT.append(ADD)
        VEC.append(D[i])
        
    MAT.append([0 for i in range(N+2)])
    MAT[-1][-1] = 1
    VEC.append(YbTrue)
    
    sweep(MAT, VEC)
    
    ys = np.array(VEC)
    ys_der = differentiate_2or_1der(ys, h)
    
    return [xs], [[np.array([ys[i], ys_der[i]]) for i in range(len(xs))]], ys_der[0]

def createGraph(method, a, b, YaTrue, YbTrue, h, Y, Y_, Name, trueSolution = None):
    XS, YS, Ya_ = method(a, b, YaTrue, YbTrue, h, Y, Y_)
    
    N = 18
    K = 9
    
    col = iter(cm.rainbow(np.linspace(0, 1, N - K)))
    
    fig, ax = plt.subplots()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    for i in range(len(XS)):
        if math.fabs(YS[i][0][1] - Ya_) < h:
            ax.plot(XS[i], [Y[0] for Y in YS[i]], '*', color = "red", label = "solution")
            break
        else:
            if i < N-1 and i>K:
                ax.plot(XS[i], [Y[0] for Y in YS[i]], color = next(col), label = "iteration " + str((i + 1-K)))
    
    if trueSolution != None:
        ax.plot(XS[i], [trueSolution(X) for X in XS[0]], '--', color = "black", label = "true solution") 
        
    ax.scatter([a, b], [YaTrue, YbTrue], label = "true", color = "green")

    ax.set_title(Name)
    ax.legend()
    plt.show()

def createGraphLinear(method, a, b, YaTrue, YbTrue, h, p, q, f, Name, trueSolution = None, deriv = False):
    XS, YS, Ya_ = method(a, b, YaTrue, YbTrue, h, p, q, f)
    col = iter(cm.rainbow(np.linspace(0, 1, len(XS))))
    
    fig, ax = plt.subplots()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    for i in range(len(XS)):
        if YS[i][0][1] == Ya_:
            ax.plot(XS[i], [Y[0] for Y in YS[i]], '.', color = "red", label = "solution")
            break
        else:
            ax.plot(XS[i], [Y[0] for Y in YS[i]], color = next(col), label = "iteration " + str((i + 1)))
    
    if trueSolution != None:
        ax.plot(XS[0], [trueSolution(X) for X in XS[0]], '--', color = "black", label = "true solution") 
        
    if deriv:
        ax.arrow(XS[0][0], YS[0][0][0], 0.2, (0.2) * YaTrue, width = 0.03, head_length = 0.04, alpha = 0.4, color = "red")
        ax.scatter([b], [YbTrue], label = "true", color = "green")
    else:
        ax.scatter([a, b], [YaTrue, YbTrue], label = "true", color = "green")

    ax.set_title(Name)
    ax.legend()
    plt.show()
    
def inflGen(method, a, b, YaTrue, YbTrue, h, Y_, Y__, trueSolution, trueSolutionDer):
    XS, YS, Ya_ = method(a, b, YaTrue, YbTrue, h, Y_, Y__)
    
    total1 = 0
    total2 = 0
    
    use = -1
    
    for i in range(len(XS)):
        if YS[i][0][1] == Ya_:
            use = i
            break
    
    if use == -1:
        print("err")
    
    for i in range(len(XS[0])):
        d1 = trueSolution(XS[0][i]) - YS[use][i][0]
        d2 = trueSolutionDer(XS[0][i]) - YS[use][i][1]
        total1 = max(total1, math.fabs(d1))
        total2 = max(total2, math.fabs(d2))
        
    return total1
    
def inflLin(method, a, b, YaTrue, YbTrue, h, p, q, f, trueSolution, trueSolutionDer):
    XS, YS, Ya_ = method(a, b, YaTrue, YbTrue, h, p, q, f)
    
    total1 = 0
    total2 = 0
    
    use = -1
    
    for i in range(len(XS)):
        if YS[i][0][1] == Ya_:
            use = i
            break
    
    if use == -1:
        print("err")
    
    for i in range(len(XS[0])):
        d1 = trueSolution(XS[0][i]) - YS[use][i][0]
        d2 = trueSolutionDer(XS[0][i]) - YS[use][i][1]
        total1 = max(total1, math.fabs(d1))
        total2 = max(total2, math.fabs(d2))
        
    return total1
    
def infls(a, b, YaTrue, YbTrue, p, q, f, Y_, Y__, trueSolution, trueSolutionDer):
    fig, ax = plt.subplots()
    
    hs = np.arange(0.25, 0.5, 0.001)
    
    inflShootingGeneral = ([inflGen(shootingMethodGeneral, a, b, YaTrue, YbTrue, h, Y_, Y__, trueSolution, trueSolutionDer) for h in hs])
    inflshootingMethodLinear = ([inflLin(shootingMethodLinear, a, b, YaTrue, YbTrue, h, p, q, f, trueSolution, trueSolutionDer) for h in hs])
    inflgridMethodLinear = ([inflLin(gridMethodLinear, a, b, YaTrue, YbTrue, h, p, q, f, trueSolution, trueSolutionDer) for h in hs])
    inflFixed = ([inflLin(fixedPointsMethodLinear, a, b, trueSolutionDer(a), YbTrue, h, p, q, f, trueSolution, trueSolutionDer) for h in hs])
    
    hs = [math.log(h) for h in hs]
    inflShootingGeneral = [math.log(j) for j in inflShootingGeneral]
    inflshootingMethodLinear = [math.log(j) for j in inflshootingMethodLinear]
    inflgridMethodLinear = [math.log(j) for j in inflgridMethodLinear]
    inflFixed = [math.log(j) for j in inflFixed]
    
    def parse(AR):
        xs = []
        ys = []
        
        for i in range(1, len(hs) - 1):
            if (AR[i - 1] < AR[i] and AR[i] > AR[i + 1]):
                xs.append(hs[i])
                ys.append(AR[i])
        return xs, ys#hs, AR
    
    def parse1(AR):
        xs = []
        ys = []
        
        for i in range(1, len(hs) - 1):
            if (AR[i - 1] > AR[i] and AR[i] < AR[i + 1]):
                xs.append(hs[i])
                ys.append(AR[i])
        return xs, ys#hs, AR
    
    xs, ys = parse1(inflShootingGeneral)
    coef = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
    print("shooting general order: " + str(round(coef)))
    ax.plot(xs, ys, label = "shooting general" + ". order: " + str(round(coef)))
    
    xs, ys = parse1(inflshootingMethodLinear)
    coef = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
    print("shooting linear order: " + str(round(coef)))
    ax.plot(xs, ys, '--', label = "shooting linear" + ". order: " + str(round(coef)))
    
    xs, ys = parse(inflgridMethodLinear)
    coef = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
    print("grid order: " + str(round(coef)))
    ax.plot(xs, ys, label = "grid" + ". order: " + str(round(coef)))
    
    xs, ys = parse(inflFixed)
    coef = (ys[-1] - ys[0]) / (xs[-1] - xs[0])
    print("fixed point order: " + str(round(coef)))
    ax.plot(xs, ys, label = "fixed points" + ". order: " + str(round(coef)))
    
    ax.set_xlabel("ln(h)")
    ax.set_ylabel("ln(Delta)")
    ax.legend()
    plt.show()




def Y_(x, y, y_):
    return y_
    
def Y__(x, y, y_):
    return -2*y_/x - y + 1/x

def p(x):
    return 2/x

def q(x):
    return 1

def f(x):
    return 1/x

def trueSol(x):
    return (-1.5*math.pi + 1) * math.sin(x)/x + math.cos(x)/x + 1/x

def trueSolDer(x):
    return (-1.5*math.pi + 1) * ( (math.sin(x) - math.cos(x)*x) / x**2 ) + (math.cos(x) + math.sin(x)*x)/x**2 - 1/x**2


a = math.pi
b = 1.5*math.pi
h = 0.001

YaTrue = 0
YbTrue = 1







createGraph(shootingMethodGeneral, a, b, YaTrue, YbTrue, h, Y_, Y__, "Shooting general", trueSol)
createGraphLinear(shootingMethodLinear, a, b, YaTrue, YbTrue, h, p, q, f, "Shooting linear method", trueSol)
createGraphLinear(gridMethodLinear, a, b, YaTrue, YbTrue, h, p, q, f, "Grid method", trueSol)
createGraphLinear(fixedPointsMethodLinear, a, b, YaTrue, YbTrue, h, p, q, f, "Fictitious point method", None, True)



def tmp():
    def Y_(x, y, y_):
        return y_
    
    def Y__(x, y, y_):
        return -2*y_/x - y + 1/x

    def p(x):
        return 2/x

    def q(x):
        return 1

    def f(x):
        return 1/x

    def trueSol(x):
        return (-1.5*math.pi + 1) * math.sin(x)/x + math.cos(x)/x + 1/x

    def trueSolDer(x):
        return (-1.5*math.pi + 1) * ( (math.sin(x) - math.cos(x)*x) / x**2 ) + (math.cos(x) + math.sin(x)*x)/x**2 - 1/x**2


    a = math.pi
    b = 1.5*math.pi
    h = 0.001

    YaTrue = 0
    YbTrue = 1

    infls(a, b, YaTrue, YbTrue, p, q, f, Y_, Y__, trueSol, trueSolDer)


tmp()