import unittest
import numpy as np
from matplotlib import pyplot as plt

import methods.Trust;
import methods.BFGS;
import methods.NM;

N = 3

def mod(v):
    return np.sqrt(np.dot(v, v))


A = np.random.uniform(size=(N, N))
b = np.random.uniform(size=N)

def f(vec):
    x = np.array(vec)
    return np.dot(A @ x - b, A @ x - b)


def fprime(vec):
    x = np.array(vec)
    return A.T @ (A @ x - b)


def hes(vec):
    x = np.array(vec)
    return A.T @ A


def NM_wrapper(epsilon, max_iters):
    starts = [np.random.uniform(size=N) for _ in range(N+1)]

    triags = methods.NM.solve(f, starts, maxiter=max_iters, epsilon=epsilon)

    return triags


def BFGS_wrapper(epsilon, max_iters):
    ps = methods.BFGS.solve(f, fprime, np.array(np.random.uniform(size=N)), maxiter=max_iters, epsi=epsilon)
    return ps


def Trust_wrapper(epsilon, max_iters):
    ps = methods.Trust.solve(f, fprime, hes, np.random.uniform(size=N), maxiter=max_iters, gtol=epsilon)
    return ps

class SinusTest(unittest.TestCase):

    def test_all(self):
        eps = 0.0001
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

        print(r1)
        print(r2)
        print(r3)

        eps = 0.01
        self.assertTrue(np.dot(r1 - r2, r1 - r2) < eps**2)
        self.assertTrue(np.dot(r3 - r2, r3 - r2) < eps**2)
        self.assertTrue(np.dot(r1 - r3, r1 - r3) < eps**2)

    def test_error(self):
        errs = np.linspace(0.00001, 0.1, 100)
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