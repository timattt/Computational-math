import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import math
import random
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from mayavi import mlab

def identity(n):
    res = []
    for i in range(n):
        tmp = []
        for j in range(n):
            tmp.append(0)
        res.append(tmp)
    return res

def Area(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

class EqPoint:
    def __init__(self):
        self.neighbors = []#[EqPoint, path_length, edge_length, normal]
        self.x = 0
        self.y = 0
        self.nearest_vor_points = []#
        self.index = 0
        self.isBorder = False
        self.area = 0
        
    def __repr__(self):
        return "point: " + str(self.index)

def createBorder(N):
    result = []
    r = 2
    for i in range(N):
        x = r*math.cos(2*math.pi * i / N)
        y = r*math.sin(2*math.pi * i / N)
        result.append([x, y])
    return result

def createRandom(N):
    result = []
    r = 2
    while len(result) != N:
        x = r * random.uniform(-1.0, 1.0)
        y = r * random.uniform(-1.0, 1.0)
        
        if (x**2 + y**2 < r*r):
            result.append([x, y])
    return result

def createBorder1(N):
    result = []
    r = 2
    for i in range(N):
        x = r*math.cos(2*math.pi * i / N)
        y = r*math.sin(2*math.pi * i / N)
        result.append([x, y])
    return result

def createRandom1(N):
    result = []
    r = 2
    for i in range(0, N):
        for j in range(0, N):
            x = -r + 2*r*i/N
            y = -r + 2*r*j/N
            if x == 0 and y == 0:
                x = 0.001
            if (x**2 + y**2 < r**2):
                result.append([x, y])
    return result

def createEqPoints(startPoints, morePoints, draw=True):
    points = startPoints + morePoints
    vor = Voronoi(points)
    
    pts = [EqPoint() for _ in range(len(vor.points))]
    
    i = 0
    for pt in vor.points:
        pts[i].x = pt[0]
        pts[i].y = pt[1]
        pts[i].index = i
        #plt.text(pt[0], pt[1], str(i))
        i += 1
    
    i = 0
    for ind in vor.ridge_points:
        start_ind = ind[0]
        end_ind = ind[1]
        
        start = pts[start_ind]
        end = pts[end_ind]
        
        ridgeStartInd = vor.ridge_vertices[i][0]
        ridgeEndInd = vor.ridge_vertices[i][1]
        
        edgeLen = 0
        if ridgeStartInd != -1 and ridgeEndInd != -1:
            edgeLen = math.sqrt((vor.vertices[ridgeStartInd][0] - vor.vertices[ridgeEndInd][0])**2 +
                                (vor.vertices[ridgeStartInd][1] - vor.vertices[ridgeEndInd][1])**2)
        else:
            start.isBorder = end.isBorder = True
        
        length = math.sqrt((start.x - end.x)**2 + (start.y - end.y)**2)
        
        norm = np.array([end.x - start.x, end.y - start.y])/length

        start.neighbors.append([end, length, edgeLen, norm])
        end.neighbors.append([start, length, edgeLen, -norm])
        
        i += 1
        
    for i in range(len(pts)):
        inds = vor.regions[vor.point_region[i]]
        if -1 in inds:
            continue
        for index in inds:
            pts[i].nearest_vor_points.append(vor.vertices[index])
    
    for pt in pts:
        pt.area = Area(pt.nearest_vor_points)
        
    if not draw:
        return pts
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot([startPoints[i][0] for i in range(-1, len(startPoints))], [startPoints[i][1] for i in range(-1, len(startPoints))])
    ax.scatter([startPoints[i][0] for i in range(-1, len(startPoints))], [startPoints[i][1] for i in range(-1, len(startPoints))], color="red")
    
    ax.scatter([morePoints[i][0] for i in range(-1, len(morePoints))], [morePoints[i][1] for i in range(-1, len(morePoints))], color="green")

    fig = voronoi_plot_2d(vor)

    for pt in pts:
        plt.fill([pt.nearest_vor_points[i][0] for i in range(len(pt.nearest_vor_points))], [pt.nearest_vor_points[i][1] for i in range(len(pt.nearest_vor_points))], color = (random.random(), random.random(), random.random()))
  
    plt.show()
        
    return pts
       
def solve(initFunc, pts, k, Vx, Vy, rightFunc, draw=True):
    N = len(pts)
    MATRIX = identity(N)
    B = []
    
    currentLine = 0
    def insertEquation(pairs, currentLine, right):
        for pair in pairs:
            index, coef = pair
            MATRIX[currentLine][index] = coef
        B.append(right)
    
    for pt in pts:
        if pt.isBorder:
            x = pt.x
            y = pt.y
            val = initFunc(x, y)
            insertEquation([[pt.index, 1]], currentLine, val)
        else:
            coef0 = 0
            x = pt.x
            y = pt.y
            for neig in pt.neighbors:
                coef0 += ((Vx * neig[3][0] + Vy * neig[3][1]) / 2 + k / neig[1]) * neig[2]
            coefs = [[pt.index, coef0]]
            for neig in pt.neighbors:
                coefs.append([neig[0].index, ((Vx * neig[3][0] + Vy * neig[3][1]) / 2 - k / neig[1]) * neig[2]])
            insertEquation(coefs, currentLine, pt.area * rightFunc(x, y))
        currentLine += 1
        
    A = np.array(MATRIX)
    b = np.array(B)
    
    U = np.linalg.solve(A, b)
    return np.array([pt.x for pt in pts]), np.array([pt.y for pt in pts]), np.array(U)
    
def makeGraph(initFunc, rightFunc, trueSol, rnd, N):
    xs, ys, Us = solve(initFunc, createEqPoints(createBorder(N), rnd(N)), -1, 0, 0, rightFunc)
    UsTrue = [trueSol(xs[i], ys[i])+4 for i in range(len(xs))]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_trisurf(xs, ys, Us, cmap = "plasma", label = "comp")
    fig.colorbar(surf)
    #surf = ax.plot_trisurf(xs, ys, UsTrue, cmap = "inferno", label = "comp")
    #fig.colorbar(surf)
    fig.tight_layout()

    plt.show()
    
def calcInfl(initFunc, rightFunc, trueSol, N):
    xs, ys, Us = solve(initFunc, createEqPoints(createBorder(N), createRandom1(N), draw=False), 1, 0, 0, rightFunc)
    UsTrue = [trueSol(xs[i], ys[i]) for i in range(len(xs))]
    
    err = 0
    for i in range(len(Us)):
        err = max(err, math.fabs(UsTrue[i] - Us[i]))
        
    return err

def makeInflGraph(initFunc, rightFunc, trueSol):
    N = np.arange(3, 15, 1)
    errs = []
    for n in N:
        errs.append(math.log(calcInfl(initFunc, rightFunc, trueSol, n)))
        
    N = np.log(N)
    
    fig = plt.figure()
    coefs = np.polyfit(N, errs, 1)
    plt.plot(N, [N[i] * coefs[0] + coefs[1] for i in range(len(N))])
    #plt.plot(N, errs)
    print(coefs)
    
    
    plt.show()

def initFunc(x, y):
    s = y / math.sqrt(x**2 + y**2)
    c = x / math.sqrt(x**2 + y**2)
    
    return s + 4/3*c+2/3*(c**2-s**2)+5

def rightFunc(x, y):
    return y**2

def trueSol(x, y):
    r = math.sqrt(x**2 + y**2)
    s = y / r
    c = x / r
    return -r**4/32 + r**4/24*(c**2-s**2) + 5.5 + r/2*s + 2*r/3*c
        
makeGraph(initFunc, rightFunc, trueSol, createRandom, 100)
makeGraph(initFunc, rightFunc, trueSol, createRandom1, 10)
makeInflGraph(initFunc, rightFunc, trueSol)