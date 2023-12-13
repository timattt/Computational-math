import methods.NM
import numpy as np
import matplotlib.pyplot as plt
import methods.BFGS
import methods.Trust
import unittest

def f(vec):
    x = vec[0]
    y = vec[1]
    return -np.exp(-x**2 - y**2)


def fprime(vec):
    x = vec[0]
    y = vec[1]
    return np.array([-2*x*f(vec), -2*y*f(vec)])


def hes(vec):
    x = vec[0]
    y = vec[1]
    fu = f(vec)
    return np.array([[-2 * fu + 4 * x**2 * fu, 4*x*y*fu],
                     [4*x*y*fu, -2*fu+4*y**2 * fu]])

def drawFunc(f, fromX, fromY, toX, toY):
    # Grid for plotting
    N = 100
    x = np.linspace(fromX, toX, num=N)
    y = np.linspace(fromY, toY, num=N)
    X, Y = np.meshgrid(x, y)

    Z = f(X, Y)

    plt.imshow(Z, extent=[fromX, toX, fromY, toY], origin="lower")
    plt.colorbar()

def drawPoints(ps):
    i = 0
    plt.plot([p[0] for p in ps], [p[1] for p in ps], color="black")
    for p in ps:
        plt.scatter([p[0]], [p[1]], label="{}".format(i), color="red")
        i += 1

def drawNM():
    v1 = np.array([0, 0])
    v2 = np.array([1, 0])
    v3 = np.array([0, 1])

    triags = methods.NM.solve(f, v1, v2, v3)

    for triag in triags:
        plt.plot([triag[0][0], triag[1][0], triag[2][0], triag[0][0]],
                 [triag[0][1], triag[1][1], triag[2][1], triag[0][1]], '--')

    plt.scatter([t[0][0] for t in triags], [t[0][1] for t in triags], color="black")
    plt.plot([t[0][0] for t in triags], [t[0][1] for t in triags], color="black")

    plt.show()

def drawTrust():
    ps = methods.Trust.solve(f, fprime, hes, [0.3, 0.4])
    drawPoints(ps)
    plt.show()
def drawBFGS():
    ps = methods.BFGS.solve(f, fprime, np.array([1, 1]))
    drawPoints(ps)
    plt.show()

class DrawTests(unittest.TestCase):

    def testDrawNM(self):
        drawFunc(lambda x, y: -np.exp(-x**2 - y**2), -0.25, -0.25, 1, 1)
        drawNM()

    def testDrawBFGS(self):
        drawFunc(lambda x, y: -np.exp(-x ** 2 - y ** 2), -0.1, -0.1, 0.05, 0.05)
        drawBFGS()

    def testDrawTrust(self):
        drawFunc(lambda x, y: -np.exp(-x ** 2 - y ** 2), -0.1, -0.1, 0.5, 0.5)
        drawTrust()





