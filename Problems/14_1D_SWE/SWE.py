import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm

g = 9.8
N = 100 # x
M = 1000 # t
h = 0.01 # dx
tau = 0.001 # dt

L = N*h

print("Task: L = {}".format(L))

H = np.zeros((N, M))
U = np.zeros((N, M))
# U[i - x][j - t]

Ustartcond = 0 * np.arange(N)
#Hstartcond = np.array([1 + np.exp(-500*(x-L/2)**2) for x in np.arange(0, L, h)])
Hstartcond = (1+h*np.exp(-(np.arange(N)-N/2)**2/10))

U[:, 0] = Ustartcond
H[:, 0] = Hstartcond

def lambda1(u, h):
    return u + np.sqrt(g * h)

def lambda2(u, h):
    return u - np.sqrt(g * h)

def lambda_plus(lam):
    return 0.5 * (lam + np.abs(lam))

def lambda_minus(lam):
    return 0.5 * (lam - np.abs(lam))

def L(F, i, j, lamb):
    return lambda_plus(lamb(U[i, j], H[i, j])) * (F[i, j] - F[i - 1, j]) / h + lambda_minus(lamb(U[i, j], H[i, j])) * (F[i+1, j] - F[i, j]) / h

def getU(i, j):
    if i < 0 or i >= N:
        return 0
    return U[i][j]

def getH(i, j):
    if i < 0 or i >= N:
        return 0
    return H[i][j]

def getC(i, j):
    if i < 0 or i >= N:
        return 0
    if H[i][j] <= 0:
        print("Bad H")
        quit(0)
    return np.sqrt(g * H[i][j])

xs = np.arange(0, N * h, h)

Hmin = 10000.0
Hmax = -10000.0

for j in range(0, M-1):
    # internal
    for i in range(1, N-1):#N - 1):
        #print("step: i={}, j={}".format(i, j))
        #print("VEC: h={}, u={}".format(H[i, j], U[i, j]))
        lam1 = lambda1(U[i, j], H[i, j])
        lam2 = lambda2(U[i, j], H[i, j])
        c = getC(i, j)
        
        #print("lambda1={}, lambda2={}, c={}".format( lam1, lam2, c))
        
        (H[i, j+1], U[i, j+1]) = tuple(np.linalg.solve(np.array([[-c/tau, H[i, j]/tau], [c/tau, H[i, j]/tau]]),
                 np.array([-c*H[i, j]/tau + c * L(H, i, j, lambda1) + H[i, j]*U[i, j]/tau - H[i, j]*L(U, i, j, lambda1),
                           c*H[i, j]/tau - c * L(H, i, j, lambda2) + H[i, j]*U[i, j]/tau - H[i, j]*L(U, i, j, lambda2)])))
        
        Hmin = min(Hmin, H[i, j+1])
        Hmax = max(Hmax, H[i, j+1])
        
        if np.abs(U[i, j+1]) > c:
            print("1. U is greater then c, i={}, j={}".format(i, j))
            print(U[i, j+1])
            print(c)
            quit(0) 
        #print("result: H={}, U={}".format(H[i, j+1], U[i, j+1]))
        
    # border
    U[0, j+1] = 0#np.sin(j / M)#U[1, j+1]
    U[-1, j+1] = 0#U[-2, j+1]
    
    # left
    lam1 = lambda1(U[0, j], H[0, j])
    lam2 = lambda2(U[0, j], H[0, j])
    c = getC(0, j)
    
    if np.abs(U[0, j]) > c:
        print("2. U is greater then c, i={}, j={}".format(i, j))
        quit(0) 
    
    #print("left border lambdas: lambda1={}, lambda2={}".format(lam1, lam2))
    
    if lam1 >= 0 and lam2 >= 0:
        print("left bad lambdas. i={}, j={}".format(0, j))
        quit(0)
        
    if lam1 < 0:
        # A * H[0, j+1] = B
        A = -c/tau
        B = -U[0, j+1]*H[0, j]/tau + -c*H[0, j]/tau + c * L(H, 0, j, lambda1) + H[0, j]*U[0, j]/tau - H[0, j]*L(U, 0, j, lambda1)
        
        if np.abs(A) < 0.0001:
            print("problem left lam1")
            quit(0)
        
        H[0, j+1] = B/A
    elif lam2 < 0:
        A = c/tau
        B = -U[0, j+1]*H[0, j] + c*H[0, j]/tau - c * L(H, 0, j, lambda2) + H[0, j]*U[0, j]/tau - H[0, j]*L(U, 0, j, lambda2)
            
        if np.abs(A) < 0.0001:
            print("problem left lam2")
            quit(0)
            
        H[0, j+1] = B/A
        
    #print("left border: U={}, H={}".format( U[0, j+1], H[0, j+1]))
    
    #print(H[:, j])
    
    #right
    lam1 = lambda1(U[-1, j], H[-1, j])
    lam2 = lambda2(U[-1, j], H[-1, j])
    c = getC(N-1, j)
    
    #print("right border lambdas: lambda1={}, lambda2={}".format(lam1, lam2))
    
    if lam1 <= 0 and lam2 <= 0:
        print("right bad lambdas. i={}, j={}".format(N-1, j))
        quit()
        
    A = 0
    B = 0
        
    if lam1 > 0:
        # A * H[-1, j] = B
        A = -c/tau
        B = -U[-1, j+1]*H[-1, j]/tau + -c*H[-1, j]/tau + c * L(H, -1, j, lambda1) + H[-1, j]*U[-1, j]/tau - H[-1, j]*L(U, -1, j, lambda1)
    
        if np.abs(A) < 0.0001:
            print("problem right lam1")
            quit(0)
    
        H[-1, j+1] = B/A
    elif lam2 < 0:
        A = c/tau
        B = -U[-1, j+1]*H[-1, j] + c*H[-1, j]/tau - c * L(H, -1, j, lambda2) + H[-1, j]*U[-1, j]/tau - H[-1, j]*L(U, -1, j, lambda2)
        
        if np.abs(A) < 0.0001:
            print("problem right lam2")
            quit(0)
            
        H[-1, j+1] = B/A
        
    #print("right border: U={}, H={}".format( U[-1, j+1], H[-1, j+1]))
    
figure, axis = plt.subplots(3)
    
axis[0].set_ylim(bottom = 0.999*np.min(H[:, 0])-h, top = 1.001*np.max(H[:, 0])+h)
axis[0].plot(xs, H[:, 0])
axis[0].set_ylabel("H(x)")
axis[0].grid()
axis[0].set_title("Start conditions")
    
axis[1].set_ylabel("H(x)")
axis[1].grid()
axis[1].set_title("Process")
axis[1].set_ylim(bottom = Hmin-h, top = Hmax+h)
line, = axis[1].plot(xs, H[:, 0])

def animate(i):
    line.set_ydata(H[:, i % M])  # update the data.
    return line,

ani = animation.FuncAnimation(figure, animate, interval=1, blit=True)

axis[2].set_title("Rainbow")
axis[2].set_xlabel("x")
axis[2].set_ylabel("H(x)")
axis[2].set_ylim(bottom = Hmin-h, top = Hmax+h)
#axis[2].set_ylim(bottom = 0.07, top = 0.13)
axis[2].set_xlim(left = 0, right = 1)
color = iter(cm.rainbow(np.linspace(0, 1, M)))

for i in range(M):
    x, y = xs, H[:, i]
    c = next(color)
    axis[2].plot(x, y, color = c)
    
plt.show()
    