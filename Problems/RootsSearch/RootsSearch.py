import math

def func(x):
    return (x-1)*(x-1)*(x-1)#math.cos(x) - x*x*x

def func_der(x):
    return 3*(x-1)*(x-1)#-math.sin(x)-3*x*x

def bin_search(a, b, epsilon, f):
    while b - a > epsilon:
        c = (a + b) / 2
        if f(b) * f(c) < 0:
            a = c
        else:
            b = c
    return (a + b) / 2

def newton_tangents(x0, epsilon, f, f_der):
    xn = x0
    xnm1 = x0-1
    q = -1
    
    for i in range(1000000):
        xnp1 = xn - f(xn) / f_der(xn)
        q = (xnp1 - xn) / (xn - xnm1)
        if math.fabs(xnp1 - xn) < epsilon:
            break
        xnm1 = xn
        xn = xnp1
    
    q = 1 / (1 - q)

    return [xn, q]

def newton_chords(x0, x1, epsilon, f):
    xn = x0
    xnm1 = x1
    q = -1
    
    for i in range(1000000):
        xnp1 = xn - f(xn) * (xn - xnm1) / (f(xn) - f(xnm1))
        q = (xnp1 - xn) / (xn - xnm1)
        if math.fabs(xnp1 - xn) < epsilon:
            break
        xnm1 = xn
        xn = xnp1
        
    q = 1 / (1 - q)
    
    return [xn, q]

def newton_tangent(x0, epsilon, f, f_der):
    xn = x0
    xnm1 = x0 - 1
    der = f_der(x0)
    q = -1
    
    for i in range(1000000):
        xnp1 = xn - f(xn) / der
        q = (xnp1 - xn) / (xn - xnm1)
        if math.fabs(xnp1 - xn) < epsilon:
            break
        xnm1 = xn
        xn = xnp1
        
    q = 1 / (1 - q)
    
    return [xn, q]

# Tests
print(bin_search(0, 1, 0.001, func))
print(newton_tangents(0.5, 0.001, func, func_der))
print(newton_tangent(0.5, 0.001, func, func_der))
print(newton_chords(200, 30, 0.001, func))

