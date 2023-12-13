import unittest
import numpy as np
from matplotlib import pyplot as plt

import methods.Trust;
import methods.BFGS;
import methods.NM;

alpha = 10 * np.sqrt(2)


def mod(v):
    return v[0] ** 2 + v[1] ** 2

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


def NM_wrapper(epsilon, max_iters):
    v1 = np.array([0, 0])
    v2 = np.array([1, 0])
    v3 = np.array([0, 1])

    triags = methods.NM.solve(f, v1, v2, v3, maxiter=max_iters, epsilon=epsilon)

    return [[t[0][0], t[0][1]] for t in triags]


def BFGS_wrapper(epsilon, max_iters):
    ps = methods.BFGS.solve(f, fprime, np.array([1, 1]), maxiter=max_iters, epsi=epsilon)
    return [[p[0], p[1]] for p in ps]


def Trust_wrapper(epsilon, max_iters):
    ps = methods.Trust.solve(f, fprime, hes, [0.3, 0.4], maxiter=max_iters, gtol=epsilon)
    return [[p[0], p[1]] for p in ps]

class SinusTest(unittest.TestCase):

    def test_all(self):
        eps = 0.00000001
        maxi = 100
        ps1 = NM_wrapper(eps, maxi)
        ps2 = BFGS_wrapper(eps, maxi)
        ps3 = Trust_wrapper(eps, maxi)

        r1 = np.array(ps1[-1])
        r2 = np.array(ps2[-1])
        r3 = np.array(ps3[-1])

        print(len(ps1))
        print(len(ps2))
        print(len(ps3))

        self.assertTrue(mod(r1-r2)**2 + mod(r1-r3)**2 + mod(r2-r3)**2 < eps ** 2)

    def test_error(self):
        errs = np.linspace(0.0000001, 0.1, 100)
        maxi = 100000

        ns1 = [len(NM_wrapper(err, maxi)) for err in errs]
        ns2 = [len(BFGS_wrapper(err, maxi)) for err in errs]
        ns3 = [len(Trust_wrapper(err, maxi)) for err in errs]

        plt.plot(errs, ns1, label="NM")
        plt.plot(errs, ns2, label="BFGS")
        plt.plot(errs, ns3, label="Trust")

        plt.xlabel("epsilon")
        plt.ylabel("iterations")
        plt.grid()
        plt.legend()
        plt.show()