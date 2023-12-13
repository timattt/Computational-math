import numpy as np


def solve(f, v1, v2, v3, alpha=1, beta=0.5, gamma=2, maxiter=5, epsilon=0.01):
    triags = []
    for i in range(maxiter):
        adict = [[v1, f(v1)], [v2, f(v2)], [v3, f(v3)]]
        points = sorted(adict, key=lambda x: x[1])

        b = points[0][0]
        g = points[1][0]
        w = points[2][0]

        mid = (g + b) / 2

        # reflection
        xr = mid + alpha * (mid - w)
        if f(xr) < f(g):
            w = xr
        else:
            if f(xr) < f(w):
                w = xr
            c = (w + mid) / 2
            if f(c) < f(w):
                w = c
        if f(xr) < f(b):

            # expansion
            xe = mid + gamma * (xr - mid)
            if f(xe) < f(xr):
                w = xe
            else:
                w = xr
        if f(xr) > f(g):

            # contraction
            xc = mid + beta * (w - mid)
            if f(xc) < f(w):
                w = xc

        # update points
        v1 = w
        v2 = g
        v3 = b

        if (v2[0] - v3[0]) ** 2 + (v2[1] - v3[1]) ** 2 < epsilon ** 2:
            break

        triags.append((v1, v2, v3))
    return triags
