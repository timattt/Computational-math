import numpy as np
import numpy.linalg as ln
import scipy as sp
import scipy.optimize

def solve(f, fprime, x0, maxiter, epsi):

    k = 0
    gfk = fprime(x0)
    N = len(x0)
    # Set the Identity matrix I.
    I = np.eye(N, dtype=int)
    Hk = I
    xk = x0

    vecs = []

    while ln.norm(gfk) > epsi and k < maxiter:
        # pk - direction of search

        pk = -np.dot(Hk, gfk)

        # Line search constants for the Wolfe conditions.
        # Repeating the line search

        # line_search returns not only alpha
        # but only this value is interesting for us

        line_search = sp.optimize.line_search(f, fprime, xk, pk)
        alpha_k = line_search[0]

        xkp1 = xk + alpha_k * pk
        sk = xkp1 - xk
        xk = xkp1

        gfkp1 = fprime(xkp1)
        yk = gfkp1 - gfk
        gfk = gfkp1

        k += 1

        ro = 1.0 / (np.dot(yk, sk))
        A1 = I - ro * sk[:, np.newaxis] * yk[np.newaxis, :]
        A2 = I - ro * yk[:, np.newaxis] * sk[np.newaxis, :]
        Hk = np.dot(A1, np.dot(Hk, A2)) + (ro * sk[:, np.newaxis] *
                                           sk[np.newaxis, :])

        vecs.append(xk)

    return vecs
