import unittest

import numpy as np
import methods.Trust
import methods.BFGS
import methods.NM

alpha = np.sqrt(2)

def mod(v):
    return v[0] ** 2 + v[1] ** 2

def f(vec):
    x = vec[0]
    y = vec[1]
    return alpha * x ** 2 + x * y + y ** 2 - 6 * x - 9 * y

def fprime(vec):
    x = vec[0]
    y = vec[1]
    return np.array([alpha * 2 * x + y - 6, x + 2 * y - 9])

def hes(vec):
    x = vec[0]
    y = vec[1]
    return np.array([[2 * alpha, 1], [1, 2]])

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
    ps = methods.Trust.solve(f, fprime, hes, [1, 4], maxiter=max_iters, gtol=epsilon)
    return [[p[0], p[1]] for p in ps]

class SimpleTest(unittest.TestCase):

    def test_all(self):
        eps = 0.0000001
        maxi = 1000
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
