#========================================================================
# imports
#========================================================================

import matplotlib.pyplot as plt
import numpy as np
import math

#========================================================================
# Differentiating
#========================================================================

def differentiate_1or_1der(points:list, step:float) -> list:
    """ Дифференцируем функцию один раз с одним порядком точности. На вход - набор значений функций и шаг аргумента между этими значениями"""
    if len(points) < 2:
        raise RuntimeError("too short array of points!")
    
    res = []
    n = len(points)
    
    for i in range(n - 1):
        res.append((points[i + 1] - points[i]) / step);
    res.append(res[len(res) - 1])
    return res
    
def differentiate_2or_1der(points:list, step:float) -> list:
    """ Дифференцируем функцию один раз с двумя порядками точности. На вход - набор значений функций и шаг аргумента между этими значениями"""
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

def differentiate_1or_2der(points:list, step:float) -> list:
    """ Дифференцируем функцию два раза с двумя порядками точности. На вход - набор значений функций и шаг аргумента между этими значениями"""
    if len(points) < 3:
        raise RuntimeError("too short array of points!")
    
    res = []
    n = len(points)
    
    # from left to right
    for i in range(n - 2):
        f0 = points[i]
        f1 = points[i+1]
        f2 = points[i+2]
        der = ( f0 - 2*f1 + f2 ) / (step * step)
        res.append(der);
    
    # right border
    for i in range(n-2, n, 1):
        f0 = points[i]
        f_1 = points[i-1]
        f_2 = points[i-2]
        der = (f_2 - 2 * f_1 + f0) / (step * step)
        res.append(der)
    
    return res
    
def differentiate_2or_2der(points:list, step:float) -> list:
    """ Дифференцируем функцию два раза с двумя порядками точности. На вход - набор значений функций и шаг аргумента между этими значениями"""
    if len(points) < 4:
        raise RuntimeError("too short array of points!")
    
    res = []
    n = len(points)
    
    # from left to right
    for i in range(n - 3):
        f0 = points[i]
        f1 = points[i+1]
        f2 = points[i+2]
        f3 = points[i+3]
        der = ( 2*f0 - 5 * f1 + 4 * f2 - f3 ) / (step * step)
        res.append(der);
    
    # right border
    for i in range(n-3, n, 1):
        f0 = points[i]
        f_1 = points[i-1]
        f_2 = points[i-2]
        f_3 = points[i-3]
        der = (2 * f0 - 5*f_1 + 4*f_2 - f_3) / (step * step)
        res.append(der)
    
    return res
    
#========================================================================
# Infilicities
#========================================================================

def infilicity(derFunc, xs:list, ys:list) -> float:
    """ Находим погрешность сеточной функции относительно непрерывного исходника. Как максимум из отдельных отклонений"""
    res = 0.0
    
    i = 0
    for x in xs:
        res = max(res, abs(derFunc(x) - ys[i]))
        i += 1
        
    return res
    
def generateInfls(func, der, diff, argStart, argEnd, stepStart, stepEnd, stepStep):
    res = [[], []]
    for step in np.arange(stepStart, stepEnd, stepStep):
        xs = np.arange(argStart, argEnd, step)
        ys = generatePointsFromArr(func, xs)
        dys = diff(ys, step)
        infl = infilicity(der, xs, dys)
        res[0].append(math.log(step))
        res[1].append(math.log(infl))
    return res
    
#========================================================================
# Generating points
#========================================================================
    
def generatePointsFromRange(func, start:float, end:float, step:float) -> list:
    """ генерируем последовательность для функции по сетке. На вход - функция и данные сетки"""
    res = []
    for x in np.arange(start, end, step):
        res.append(func(x))
    return res

def generatePointsFromArr(func, xs:list) -> list:
    """ генерируем последовательность для функции по сетке. На вход - функция и сама сетка"""
    res = []
    for x in xs:
        res.append(func(x))
    return res
    
#========================================================================
# Drawing
#========================================================================
    
def drawGraph4a(xs:list, ys1:list, label1, ys2:list, label2, ys3:list, label3, ys4:list, label4):
    """ рисуем график из четырех линий"""
    fig, ax = plt.subplots()
    
    ax.plot(xs, ys1, label = label1)
    ax.plot(xs, ys2, label = label2)
    ax.plot(xs, ys3, label = label3)
    ax.plot(xs, ys4, label = label4)

    ax.set_xlabel("ln(step)")
    ax.set_ylabel("ln(infilicity)")
    
    ax.grid()
    ax.legend()
    
    plt.show()

#========================================================================
# Input
#========================================================================

stepStart = 0.0002
stepEnd = 0.01
stepStep = 0.0001

xStart = 0.2
xEnd = 1.0

def f(x:float) -> float:
    return x*x*x*x
    
def f_der(x:float) -> float:
    return 4*x*x*x
    
def f_der_der(x:float) -> float:
    return 12 * x * x;

#========================================================================
# Main part
#========================================================================

res1 = generateInfls(f, f_der, differentiate_1or_1der, xStart, xEnd, stepStart, stepEnd, stepStep)
res2 = generateInfls(f, f_der, differentiate_2or_1der, xStart, xEnd, stepStart, stepEnd, stepStep)
res3 = generateInfls(f, f_der_der, differentiate_1or_2der, xStart, xEnd, stepStart, stepEnd, stepStep)
res4 = generateInfls(f, f_der_der, differentiate_2or_2der, xStart, xEnd, stepStart, stepEnd, stepStep)
drawGraph4a(res1[0], res1[1], "1 der, 1 order", res2[1], "1 der, 2 order", res3[1], "2 der, 1 order", res4[1], "2 der, 2 order")