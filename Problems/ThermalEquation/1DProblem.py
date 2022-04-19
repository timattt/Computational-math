import math
import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

def identity(n):
    res = []
    for i in range(n):
        tmp = []
        for j in range(n):
            tmp.append(0)
        res.append(tmp)
    return res

def SwapRows(A, B, row1, row2):
    A[row1], A[row2] = A[row2], A[row1]
    B[row1], B[row2] = B[row2], B[row1]

def DivideRow(A, B, row, divider):
    A[row] = [a / divider for a in A[row]]
    B[row] /= divider

def CombineRows(A, B, row, source_row, weight):
    A[row] = [(a + k * weight) for a, k in zip(A[row], A[source_row])]
    B[row] += B[source_row] * weight

def gauss_simple(A, B):
    res = identity(len(A))
    n = len(A)
    
    for col in range(n):
        
        # Выбор строки с ненулевым элементом
        #=======================================
        if A[col][col] == 0:
            toSwap = None
            for row in range(col, n):
                if A[row][col] != 0:
                    toSwap = row
                    break
            if toSwap is None:
                return 0

            SwapRows(A, B, col, toSwap)
        #=======================================
        
        
        # Нормируем
        #=======================================
        DivideRow(A, B, col, A[col][col])
        #=======================================
        
        
        # Обрабатываем другие строки
        #=======================================
        for row in range(0, n):
            if row != col:
                CombineRows(A, B, row, col, -A[row][col])
        #=======================================
    
    return 1

def gauss_lead(A, B):
    res = identity(len(A))
    n = len(A)
    
    for col in range(n):
        
        # Выбор строки с наибольшим элементом
        #=======================================
        max_ = 0
        toSwap = None
        for row in range(col, n):
            if toSwap is None or abs(A[row][col]) > max_:
                max_ = abs(A[row][col])
                toSwap = row
        if toSwap == None:
            return 0
        #=======================================    


        # Меняем изначальную строку на нужную
        #=======================================
        SwapRows(A, B, col, toSwap)
        #=======================================
        
        
        # Нормируем
        #=======================================
        DivideRow(A, B, col, A[col][col])
        #=======================================
        
        
        # Обрабатываем другие строки
        #=======================================
        for row in range(0, n):
            if row != col:
                CombineRows(A, B, row, col, -A[row][col])
        #=======================================
        
    return 1

def FancyPrint(A, B, selected):
    for row in range(len(B)):
        print("(", end='')
        for col in range(len(A[row])):
             print("\t{1:10.2f}{0}".format(" " if (selected is None or selected != (row, col)) else "*", A[row][col]), end='')
        print("\t) * (\tX{0}) = (\t{1:10.2f})".format(row + 1, B[row]))
    print("\n")

def solve(N, sigma, h, tau, Ux0, U0t, ULt):
    alpha = sigma * tau / h**2
    beta = (1-sigma)*tau/h**2
    
    MATRIX = identity(N*N)
    B = []
    
    currentLine = 0
    
    def toMat(n, i):
        return n * N + i
    def fromMat(index):
        return (index / N, index % N)
    def insertEquation(pairs, currentLine):
        for pair in pairs:
            index, coef = pair
            MATRIX[currentLine][index] = coef
        B.append(0)
    
    for i in range(0, N):
        MATRIX[currentLine][toMat(i, 0)] = 1
        currentLine += 1
        B.append(U0t(i * tau))

    for i in range(0, N):    
        MATRIX[currentLine][toMat(i, N-1)] = 1
        currentLine += 1
        B.append(ULt(i * tau))
        
    for i in range(1, N-1):
        MATRIX[currentLine][toMat(0, i)] = 1
        currentLine += 1
        B.append(Ux0(i * h))
    
    for i in range(1, N-1):
        for n in range(0, N-1):
            insertEquation([
                    (toMat(n+1, i+1), beta),
                    (toMat(n+1, i), -1-2*beta),
                    (toMat(n+1,i-1), beta),
                    (toMat(n, i+1), alpha),
                    (toMat(n, i), 1-2*alpha),
                    (toMat(n,i-1), alpha)
                ], currentLine)
            currentLine += 1
    
    A = np.array(MATRIX)
    b = np.array(B)
    x = np.linalg.solve(A, b)
    
    U = [[0 for _ in range(N)] for _ in range(N)]

    for i in range(N):
        for n in range(N):
            U[n][i] = x[toMat(n, i)]
    
    return U

def plotrod(u_k, k, time):
    # Clear the current plot figure
    plt.clf()

    #plt.xlabel("x")
    #plt.ylabel("T(x)")

    # This is to plot u_k (u at time-step k)
    #plt.plot(k, u_k, label = "time = " + str(time))
    #plt.legend()
    
    
    plt.xlabel("x")
    plt.ylabel("y")
    tmp = [[0 for _ in range(len(u_k))], [0 for _ in range(len(u_k))], u_k, [0 for _ in range(len(u_k))], [0 for _ in range(len(u_k))]]
    plt.pcolormesh(tmp, cmap="Reds", vmin=np.min(u_k), vmax=np.max(u_k))
    plt.plot([0, max(k)], [2, 2], color = "black")
    plt.plot([0, max(k)], [3, 3], color = "black")
    
    return plt

def plotgraph(u_k, k, time):
    # Clear the current plot figure
    plt.clf()

    plt.xlabel("x")
    plt.ylabel("T(x)")

    # This is to plot u_k (u at time-step k)
    plt.plot(k, u_k)
    
    return plt

def draw(N, sigma, h, tau, Ux0, U0t, ULt):
    u = solve(N, sigma, h, tau, Ux0, U0t, ULt)
    xs = np.linspace(0, N * h, N)
    print(u)
    max_iter_time = 750
    #anim = animation.FuncAnimation(plt.figure(), lambda k: plotrod(u[k % N], xs, k), interval=1, frames=max_iter_time, repeat=True)
    #anim.save('animation1.gif', fps=25)
    #plt.show()
    anim = animation.FuncAnimation(plt.figure(), lambda k: plotgraph(u[k % N], xs, k), interval=1, frames=max_iter_time, repeat=True)
    plt.show()
    anim.save('animation2.gif', fps=25)

def Ux01(x):
    return math.exp(-x**2)
def U0t1(t):
    return 1/math.sqrt(1+4*t)
def ULt1(t):
    return math.exp(-1/(1+4*t)) / math.sqrt(1+4*t)
def uTrue(t, x):
    return 1/math.sqrt(1+4*t) * math.exp(-x**2 / (1+4*t))


def getInfl(N, sigma, h, tau):
    u = solve(N, sigma, h, tau, Ux01, U0t1, ULt1)

    res = 0
    
    for i in range(N):
        for n in range(N):
            res = max(res, math.fabs(u[n][i] - uTrue(n * tau, i * h)))
            
    return res

def calcInfls(sigma, N):
    infls = []
    if sigma > 0.5:
        sigma = 0.5
    for n in range(5, N):
        infls.append(math.log(getInfl(n, sigma, 1.0/n, 1.0/n)))
        
    if sigma == 0.5:
        return np.array(infls) * 2
    return np.array(infls)

def drawInfls(N):
    ns = [math.log(n) for n in range(5, N)]
    infls1 = calcInfls(0, N)
    infls2 = calcInfls(0.5, N)
    infls3 = calcInfls(1, N)

    fig, ax = plt.subplots()
    
    ax.plot(ns, infls1, label = "sigma = 0")
    ax.plot(ns, infls2, label = "sigma = 0.5")
    ax.plot(ns, infls3, label = "sigma = 1")

    print(round(-(infls1[-1]-infls1[0])/(ns[-1]-ns[0])))
    print(round(-(infls2[-1]-infls2[0])/(ns[-1]-ns[0])))
    print(round(-(infls3[-1]-infls3[0])/(ns[-1]-ns[0])))

    plt.legend()
    plt.show()

def Ux0(x):
    return 0
def U0t(t):
    return 10
def ULt(t):
    return 10

h = 1
tau = 2
N = 50
sigma = 0

draw(N, sigma, h, tau, Ux0, U0t, ULt)
#drawInfls(20)