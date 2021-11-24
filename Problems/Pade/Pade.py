import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.linalg import null_space

def printPol(coefs):
    items = []
    for i, x in enumerate((coefs)):
        if not x:
            continue
        items.append('{}x^{}'.format("{0:0.2f}".format(x) if x != 1 or i == 0 else '', i))
    result = ' + '.join(items)
    result = result.replace('x^0', '')
    result = result.replace('^1 ', ' ')
    result = result.replace('+ -', '- ')
    
    return result

def ran(a, b):
    if a < b:
        return range(a, b + 1)
    if a == b:
        return range(a, a + 1)
    if b < a:
        return range(b, a + 1)



def Pade(coefs, n, m):
    M = n - m + 1
    N = n + m
    T = [np.matrix([[coefs(i + j - k) for i in ran(k, M)] for j in ran(k, N)]) for k in ran(M, N)]
    r = [np.linalg.matrix_rank(T[k-M]) for k in range(M, N + 1)]
    d = [0] + [k - M + 1 - r[k-M] for k in range(M, N + 1)] + [N-M+2]
    delta = [d[k - M + 1] - d[k-M] for k in range(M, N+1)]
    
    mu1 = 0
    while delta[mu1] == 0:
        mu1 += 1
    mu2 = mu1 + 1
    while delta[mu2] == 1:
        mu2 += 1

    g = null_space(T[mu1 + 1])
    q = np.array(g[0])
    s = len(q)
    M = np.matrix([[coefs(i-j) for j in range(s)] for i in range(n+2)])
    p = M.dot(q)
    
    return np.squeeze(np.asarray(p)), q

def coefs_test(k):
    if k < 0:
        return 0
    else:
        return 1.0/math.factorial(k)
    
P, Q = Pade(coefs_test, 3, 3)

fig, ax = plt.subplots()

xs = []
ys1 = []
ys2 = []

for x in np.arange(-2, 5, 0.1):
    xs.append(x)
    ys1.append(math.exp(x))
    
    chis = 0
    po = 1
    for i in range(len(P)):
        chis += P[i] * po
        po *= x
    znam = 0
    po = 1
    for i in range(len(Q)):
        znam += Q[i] * po
        po *= x
    ys2.append(chis/znam)

chis = printPol(P)
znam = printPol(Q)
brace = ""
print(chis)
for i in range(max(len(chis), len(znam))):
    brace += "-"
print(brace)
print(znam)

plt.plot(xs, ys1, color="green", label="original\n")
plt.plot(xs, ys2, color="red", label="Pade:\n\n" + chis + "\n" + brace + brace + "\n" + znam)
ax.set_xlabel("x")
ax.set_ylabel("y")
    
ax.grid()
ax.legend()
plt.show()
