import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm

g = 9.8
N = 30 # x and y
M = 200 # t
h = 0.01 # dx and dy
tau = 0.001 # dt

L = N*h

print("Task: L = {}".format(L))

H1 = np.zeros((N, N))
H2 = np.zeros((N, N))
Q1 = np.zeros((N, N))
Q2 = np.zeros((N, N))
P1 = np.zeros((N, N))
P2 = np.zeros((N, N))

H = np.zeros((N, N, M))
U = np.zeros((N, N, M))
V = np.zeros((N, N, M))
P = np.zeros((N, N, M))
Q = np.zeros((N, N, M))
# U[i - x][j - y][k - t]

Ustartcond = 0 * np.zeros((N, N))
Vstartcond = 0 * np.zeros((N, N))

# test1
Hmin = 1
Hmax = 1 + 0.3*h
Hanimmax = 1.0005
Hstartcond = np.array([[1 + 0.3*h*np.exp(-500*((x-L/2)**2 + (y-L/2)**2)) for y in np.arange(0, L, h)] for x in np.arange(0, L, h)])

#test2
#Hmin = 1-0.3*h
#Hmax = 1 + 0.3*h
#Hanimmax = 1.0005
#Hstartcond = np.array([[1 + (h/10 if x > L/2 else 0) for y in np.arange(0, L, h)] for x in np.arange(0, L, h)])


U[:, :, 0] = Ustartcond
V[:, :, 0] = Vstartcond
H[:, :, 0] = Hstartcond

def lambdaA1(H, u, v):
    return u

def lambdaA2(H, u, v):
    return u - np.sqrt(g * H)

def lambdaA3(H, u, v):
    return u + np.sqrt(g * H)

def lambdaB1(H, u, v):
    return v

def lambdaB2(H, u, v):
    return v - np.sqrt(g * H)

def lambdaB3(H, u, v):
    return v + np.sqrt(g * H)

def lambda_plus(lam):
    return 0.5 * (lam + np.abs(lam))

def lambda_minus(lam):
    return 0.5 * (lam - np.abs(lam))

def Lx(F, i, j, k, lamb):
    lp = lambda_plus(lamb(H[i, j, k], U[i, j, k], V[i, j, k]))
    lm = lambda_minus(lamb(H[i, j, k], U[i, j, k], V[i, j, k]))
    
    res = 0
    
    if lp != 0:
        res += lp * (F[i, j, k] - F[i - 1, j, k]) / h
        
    if lm != 0:
        res += lm * (F[i+1, j, k] - F[i, j, k]) / h
    
    return res

def Ly(F, i, j, k, lamb):
    lp = lambda_plus(lamb(H[i, j, k], U[i, j, k], V[i, j, k]))
    lm = lambda_minus(lamb(H[i, j, k], U[i, j, k], V[i, j, k]))
    
    res = 0
    
    if lp != 0:
        res += lp * (F[i, j, k] - F[i, j - 1, k]) / h
        
    if lm != 0:
        res += lm * (F[i, j+1, k] - F[i, j, k]) / h
    
    return res

def getU(i, j, k):
    if i < 0 or i >= N:
        return 0
    if j < 0 or j >= N:
        return 0
    return U[i, j, k]

def getV(i, j, k):
    if i < 0 or i >= N:
        return 0
    if j < 0 or j >= N:
        return 0
    return V[i, j, k]

def getH(i, j, k):
    if i < 0 or i >= N:
        return 0
    if j < 0 or j >= N:
        return 0
    return H[i, j, k]

def getC(i, j, k):
    if i < 0 or i >= N:
        return 0
    if j < 0 or j >= N:
        return 0
    if H[i, j, k] <= 0:
        print("Bad H")
        quit(0)
    return np.sqrt(g * H[i, j, k])

for i in range(N):
    for j in range(N):
        P[i, j, 0] = V[i, j, 0] * H[i, j, 0]
        Q[i, j, 0] = U[i, j, 0] * H[i, j, 0]

for k in range(0, M-1):
    
    # A
    for j in range(N):
        #internal
        for i in range(1, N-1):
            c = getC(i, j, k)
            u = getU(i, j, k)
            v = getV(i, j, k)
            
            l1 = np.array([0, 0, 1])
            l2 = np.array([1, u-c, v])
            l3 = np.array([1, u+c, v])
            
            # Ht, qt, pt
            (Ht, Qt, Pt) = tuple(np.linalg.solve(
                np.array([[0, 0, 1],
                          [1, u-c, v],
                          [1, u+c, v]]),
                 np.array([-np.dot(l1, np.array([Lx(H, i, j, k, lambdaA1), Lx(Q, i, j, k, lambdaA1), Lx(P, i, j, k, lambdaA1)])),
                           -np.dot(l2, np.array([Lx(H, i, j, k, lambdaA2), Lx(Q, i, j, k, lambdaA2), Lx(P, i, j, k, lambdaA2)])),
                           -np.dot(l3, np.array([Lx(H, i, j, k, lambdaA3), Lx(Q, i, j, k, lambdaA3), Lx(P, i, j, k, lambdaA3)]))])))
            
            H1[i, j] = H[i, j, k] + Ht * tau
            Q1[i, j] = Q[i, j, k] + Qt * tau
            P1[i, j] = P[i, j, k] + Pt * tau
            
        # border
        Q1[0, j] = P1[0, j] = Q1[-1, j] = P1[-1, j] = 0
        
        # left
        lam2 = lambdaA2(H[0, j, k], U[0, j, k], V[0, j, k])
        lam3 = lambdaA3(H[0, j, k], U[0, j, k], V[0, j, k])
        c = np.sqrt(H[0, j, k] * g)
        u = U[0, j, k]
        v = V[0, j, k]
        
        l1 = np.array([0, 0, 1])
        l2 = np.array([1, u-c, v])
        l3 = np.array([1, u+c, v])
        
        if lam2 >= 0 and lam3 >= 0:
            print("left bad lambdas. i={}, j={}, k={}".format(0, j, k))
            quit(0)
        
        if lam2 < 0:
            Ht = -np.dot(l2, np.array([Lx(H, 0, j, k, lambdaA2), Lx(Q, 0, j, k, lambdaA2), Lx(P, 0, j, k, lambdaA2)]))
            H1[0, j] = H[0, j, k] + Ht*tau
        elif lam3 < 0:
            Ht = -np.dot(l3, np.array([Lx(H, 0, j, k, lambdaA3), Lx(Q, 0, j, k, lambdaA3), Lx(P, 0, j, k, lambdaA3)]))
            H1[0, j] = H[0, j, k] + Ht*tau  
            
        # right
        lam2 = lambdaA2(H[-1, j, k], U[-1, j, k], V[-1, j, k])
        lam3 = lambdaA3(H[-1, j, k], U[-1, j, k], V[-1, j, k])
        c = np.sqrt(H[-1, j, k] * g)
        u = U[-1, j, k]
        v = V[-1, j, k]
        
        l1 = np.array([0, 0, 1])
        l2 = np.array([1, u-c, v])
        l3 = np.array([1, u+c, v])
        
        if lam2 <= 0 and lam3 <= 0:
            print("right bad lambdas. i={}, j={}, k={}".format(-1, j, k))
            quit(0)
        
        if lam2 > 0:
            Ht = -np.dot(l2, np.array([Lx(H, -1, j, k, lambdaA2), Lx(Q, -1, j, k, lambdaA2), Lx(P, -1, j, k, lambdaA2)]))
            H1[-1, j] = H[-1, j, k] + Ht*tau
        elif lam3 > 0:
            Ht = -np.dot(l3, np.array([Lx(H, -1, j, k, lambdaA3), Lx(Q, -1, j, k, lambdaA3), Lx(P, -1, j, k, lambdaA3)]))
            H1[-1, j] = H[-1, j, k] + Ht*tau
            
    for i in range(N):
        for j in range(N):
            H[i, j, k] = H1[i, j]
            Q[i, j, k] = Q1[i, j]
            P[i, j, k] = P1[i, j]
            U[i, j, k] = Q[i, j, k] / H[i, j, k]
            V[i, j, k] = P[i, j, k] / H[i, j, k]   
            
    # B
    for i in range(N):
        #internal
        for j in range(1, N-1):
            c = getC(i, j, k)
            u = getU(i, j, k)
            v = getV(i, j, k)
            
            l1 = np.array([0, 1, 0])
            l2 = np.array([1, u, v-c])
            l3 = np.array([1, u, v+c])
            
            # Ht, qt, pt
            (Ht, Qt, Pt) = tuple(np.linalg.solve(
                np.array([[0, 1, 0],
                          [1, u, v-c],
                          [1, u, v+c]]),
                 np.array([-np.dot(l1, np.array([Ly(H, i, j, k, lambdaB1), Ly(Q, i, j, k, lambdaB1), Ly(P, i, j, k, lambdaB1)])),
                           -np.dot(l2, np.array([Ly(H, i, j, k, lambdaB2), Ly(Q, i, j, k, lambdaB2), Ly(P, i, j, k, lambdaB2)])),
                           -np.dot(l3, np.array([Ly(H, i, j, k, lambdaB3), Ly(Q, i, j, k, lambdaB3), Ly(P, i, j, k, lambdaB3)]))])))
            
            H2[i, j] = H[i, j, k] + Ht * tau
            Q2[i, j] = Q[i, j, k] + Qt * tau
            P2[i, j] = P[i, j, k] + Pt * tau
            
        # border
        Q2[i, 0] = P2[i, 0] = Q2[i, -1] = P2[i, -1] = 0
        
        # left
        lam2 = lambdaB2(H[i, 0, k], U[i, 0, k], V[i, 0, k])
        lam3 = lambdaB3(H[i, 0, k], U[i, 0, k], V[i, 0, k])
        c = np.sqrt(H[i, 0, k] * g)
        u = U[i, 0, k]
        v = V[i, 0, k]
        l1 = np.array([0, 1, 0])
        l2 = np.array([1, u, v-c])
        l3 = np.array([1, u, v+c])
        
        if lam2 >= 0 and lam3 >= 0:
            print("left bad lambdas. i={}, j={}, k={}".format(0, j, k))
            quit(0)
        
        if lam2 < 0:
            Ht = -np.dot(l2, np.array([Ly(H, i, 0, k, lambdaA2), Ly(Q, i, 0, k, lambdaA2), Ly(P, i, 0, k, lambdaA2)]))
            H2[i, 0] = H[i, 0, k] + Ht*tau
        elif lam3 < 0:
            Ht = -np.dot(l3, np.array([Ly(H, i, 0, k, lambdaA3), Ly(Q, i, 0, k, lambdaA3), Ly(P, i, 0, k, lambdaA3)]))
            H2[i, 0] = H[i, 0, k] + Ht*tau  
            
        # right
        lam2 = lambdaA2(H[i, -1, k], U[i, -1, k], V[i, -1, k])
        lam3 = lambdaA3(H[i, -1, k], U[i, -1, k], V[i, -1, k])
        c = np.sqrt(H[i, -1, k] * g)
        u = U[i, -1, k]
        v = V[i, -1, k]
        l1 = np.array([0, 1, 0])
        l2 = np.array([1, u, v-c])
        l3 = np.array([1, u, v+c])
        
        if lam2 <= 0 and lam3 <= 0:
            print("right bad lambdas. i={}, j={}, k={}".format(i, -1, k))
            quit(0)
        
        if lam2 > 0:
            Ht = -np.dot(l2, np.array([Ly(H, i, -1, k, lambdaA2), Ly(Q, i, -1, k, lambdaA2), Ly(P, i, -1, k, lambdaA2)]))
            H2[i, -1] = H[i, -1, k] + Ht*tau
        elif lam3 > 0:
            Ht = -np.dot(l3, np.array([Ly(H, i, -1, k, lambdaA3), Ly(Q, i, -1, k, lambdaA3), Ly(P, i, -1, k, lambdaA3)]))
            H2[i, -1] = H[i, -1, k] + Ht*tau  
            
    for i in range(N):
        for j in range(N):
            Q[i, j, k+1] = Q2[i, j]
            P[i, j, k+1] = P2[i, j]
            H[i, j, k+1] = H2[i, j]
            U[i, j, k+1] = Q[i, j, k+1] / H[i, j, k+1]
            V[i, j, k+1] = P[i, j, k+1] / H[i, j, k+1]
            
            if V[i, j, k+1]**2 > g * H[i, j, k+1] or U[i, j, k+1]**2 > g * H[i, j, k+1]:
                print("Bad U, V. i={}, j={}, k={}".format(i, j, k))
                quit(0)
     

# start
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = y = np.arange(0, L, h)
X, Y = np.meshgrid(x, y)
zs = Hstartcond
Z = zs.reshape(X.shape)

sur = ax.plot_surface(X, Y, Z, vmin=Hmin, vmax=Hmax, cmap=cm.coolwarm)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_zlim(Hmin, Hmax)
plt.show()
    
# process       
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = y = np.arange(0, L, h)
X, Y = np.meshgrid(x, y)
zs = H[:, :, 0]
Z = zs.reshape(X.shape)

sur = ax.plot_surface(X, Y, Z, vmin=Hmin, vmax=Hanimmax, cmap=cm.coolwarm)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_zlim(Hmin, Hmax)

def update(frame):
    ax.cla()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_zlim(1,1+h)
    zs = H[:, :, frame % M]
    Z = zs.reshape(X.shape)
    ax.plot_surface(X, Y, Z, vmin=Hmin, vmax=Hanimmax, cmap=cm.coolwarm)
    plt.title("T={}/{}".format(frame, M))

    return fig,

anim = animation.FuncAnimation(fig = fig, func = update, frames = M, interval = 4, repeat = True)

plt.show()

anim.save("test.gif")