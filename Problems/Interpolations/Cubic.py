import matplotlib.pyplot as plt
import numpy as np
import math

class Spline:
    def __init__(self, a, b, c, d, x):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.x = x
 

def BuildSpline(x, y):
    n = len(x)
    
    splines = [Spline(0, 0, 0, 0, 0) for i in range(0, n)]
    for i in range(0, n):
        splines[i].x = x[i]
        splines[i].a = y[i]
    
    splines[0].c = 0.0
    splines[n - 1].c = 0.0
    
    # SLAE
    alpha = [0.0 for i in range(0, n - 1)]
    beta = [0.0 for i in range(0, n - 1)]
 
    for i in range(1, n - 1):
        hi  = x[i] - x[i - 1]
        hi1 = x[i + 1] - x[i]
        A = hi
        C = 2.0 * (hi + hi1)
        B = hi1
        F = 6.0 * ((y[i + 1] - y[i]) / hi1 - (y[i] - y[i - 1]) / hi)
        z = (A * alpha[i - 1] + C)
        alpha[i] = -B / z
        beta[i] = (F - A * beta[i - 1]) / z
  
 
    # reverse
    for i in range(n - 2, 0, -1):
        splines[i].c = alpha[i] * splines[i + 1].c + beta[i]
    
    # find other
    for i in range(n - 1, 0, -1):
        hi = x[i] - x[i - 1]
        splines[i].d = (splines[i].c - splines[i - 1].c) / hi
        splines[i].b = hi * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 + (y[i] - y[i - 1]) / hi
        
    splines[0] = splines[1]
    return splines
 
 
def Cubic(splines, x):
    n = len(splines)
    s = Spline(0, 0, 0, 0, 0)
    
    if x <= splines[0].x: 
        s = splines[0]
    elif x >= splines[n - 1].x:
        s = splines[n - 1]
    else:
        for i in range(len(splines)):
            if splines[i].x > x:
                s = splines[i]
                break
    
    dx = x - s.x
    return s.a + (s.b + (s.c / 2.0 + s.d * dx / 6.0) * dx) * dx;
    
 
x = [1, 3, 7, 9]
y = [5, 6, 7, 8]
 
spline = BuildSpline(x, y)
xs = np.arange(-2, 10, 0.01)
fs = [Cubic(spline, x) for x in xs]
fig, ax = plt.subplots()
plt.scatter(x, y, color="green", label="origins")
plt.plot(xs, fs, label="interpolated")
ax.set_xlabel("x")
ax.set_ylabel("y")
    
ax.grid()
ax.legend()
plt.show()
