import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm

def task(g, N, M, L, T, Ustart, Hstart, Uleftbord, Urightbord):
    h = L/N
    tau = T/M
    
    H = np.zeros((N, M))
    U = np.zeros((N, M))
    # U[i - x][j - t]
    
    xs = np.arange(0, N * h, h)
    
    Ustartcond = np.array([Ustart(x) for x in xs])
    Hstartcond = np.array([Hstart(x) for x in xs])
    
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
        U[0, j+1] = Uleftbord(tau*(j+1))#np.sin(j / M)#U[1, j+1]
        U[-1, j+1] = Urightbord(tau*(j+1))#U[-2, j+1]
        H[0, j+1] = H[-1, j+1] = 0.01
        continue
        
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

    return U, H

def test(N, M):
    L = 1
    T = 1
    
    alpha = 0.000001
    g = 1
    
    def Utrue(x, t):
        return alpha*( np.sin( np.pi/L*(x-np.sqrt(g)*t) - np.sin(np.pi/L*(x+np.sqrt(g)*t)) )  )/(2*np.sqrt(g))
    def Htrue(x, t):
        return 0.01+alpha*( np.sin( np.pi/L*(x-np.sqrt(g)*t) + np.sin(np.pi/L*(x+np.sqrt(g)*t)) )  )/2
    
    def Ustart(x):
        return Utrue(x, 0)
    
    def Hstart(x):
        return Htrue(x, 0)
    
    def Uleftbord(t):
        return Utrue(0, t)

    def Urightbord(t):
        return Utrue(L, t)

    U, H = task(g, N, M, L, T, Ustart, Hstart, Uleftbord, Urightbord)

    err = 0
    for i in range(N):
        for j in range(M):
            if i > j and N - i > j:
                err += abs(Utrue(i*L/N, j*T/M) - U[i, j])
                err += abs(Htrue(i*L/N, j*T/M) - H[i, j])
    
    return err

ns = []
errs = []

for N in range(30, 100, 10):
    print(N)
    ns.append(np.log(N))
    errs.append(np.log(test(30, N)))

coefs = np.polyfit(ns, errs, 1)
print("order={}".format(-round(coefs[0])))
plt.plot(ns, [coefs[0] * n + coefs[1] for n in ns], label = "solution", color = "blue")
plt.scatter(ns, errs, label = "linear aprox", color = "red")
plt.xlabel("log(M)")
plt.ylabel("log(err)")
plt.legend()
plt.grid()
plt.show()