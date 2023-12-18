import numpy as np


def solve(f, startPoints, alpha=1, beta=0.5, gamma=2, maxiter=5, epsilon=0.01):
    triags = []
    for i in range(maxiter):
        adict = [[v, f(v)] for v in startPoints]
        points = sorted(adict, key=lambda x: x[1])

        b = points[0][0]
        g = points[1][0]
        w = points[-1][0]
        tmp = [p[0] for p in points[:-1]]
        mid = np.sum(tmp, axis=0) / (len(points) - 1)

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
        points[0][0] = b
        points[1][0] = g
        points[-1][0] = w

        startPoints = [p[0] for p in points]

        if np.dot(g - b, g - b) <= epsilon ** 4:
            break

        triags.append(b)
    return triags
