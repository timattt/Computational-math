import math
import copy
from Tools.scripts.pdeps import inverse

#=================================================================
def identity(n):
    res = []
    for i in range(n):
        tmp = []
        for j in range(n):
            if i == j:
                tmp.append(1)
            else:
                tmp.append(0)
        res.append(tmp)
    return res

def inputMatrix(n):
    return [[float(elem) for elem in input().split()] for i in range(n)]
    
def inputVector(n):
    return [float(elem) for elem in input().split()]

def FancyPrint(A, B, selected):
    for row in range(len(B)):
        print("(", end='')
        for col in range(len(A[row])):
             print("\t{1:10.2f}{0}".format(" " if (selected is None
or selected != (row, col)) else "*", A[row][col]), end='')
        print("\t) * (\tX{0}) = (\t{1:10.2f})".format(row + 1, B[row]))
    print("\n")

def SwapRows(A, B, row1, row2):
    A[row1], A[row2] = A[row2], A[row1]
    B[row1], B[row2] = B[row2], B[row1]

def DivideRow(A, B, row, divider):
    A[row] = [a / divider for a in A[row]]
    B[row] /= divider

def CombineRows(A, B, row, source_row, weight):
    A[row] = [(a + k * weight) for a, k in zip(A[row], A[source_row])]
    B[row] += B[source_row] * weight
    
def matrixmul(A, B):
    C = [[0 for row in range(len(A))] for col in range(len(B[0]))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                C[i][j] += A[i][k]*B[k][j]
    return C
    
def matrixsub(A, B):
    C = [[0 for row in range(len(A))] for col in range(len(B[0]))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            C[i][j] = A[i][j] - B[i][j]
            
    return C
    
def vectorsum(A, B):
    C = []
    
    for i in range(len(A)):
        C.append(A[i] + B[i])
        
    return C
    
def vectorsub(A, B):
    C = []
    
    for i in range(len(A)):
        C.append(A[i] - B[i])
        
    return C
    
def delta(A, B):
    res = 0
    for i in range(len(A)):
        res = max(abs(A[i] - B[i]), res)
    return res
    
def vectormul(A, B):
    C = [0 for row in range(len(A))]
    
    for i in range(len(B)):
        v = 0
        for j in range(len(B)):
            v += (B[j] * A[i][j])
        C[i] = v
        
    return C
    
def vectormulf(A, b):
    C = [0 for row in range(len(A))]
    
    for i in range(len(b)):
        C[i] = A[i] * b
        
    return C

def dot(A, B):
    res = 0
    for i in range(len(A)):
        res += A[i] * B[i]
    return res

def veccopy(A):
    C = [0 for row in range(len(A))]
    
    for i in range(len(A)):
        C[i] = A[i]
        
    return C

def matrix_from_multidimensional_function(F, X):
    n = len(F)
    res = [[0 for row in range(n)] for col in range(n)]
    
    for i in range(n):
        for j in range(n):
            res[i][j] = F[i][j](X)
            
    return res

def vector_from_multidimensional_function(F, X):
    n = len(F)
    res = [0 for row in range(n)]
    
    for i in range(n):
        res[i] = F[i](X)
            
    return res

def minor(A, i, j):
    M = copy.deepcopy(A)
    del M[i]
    for i in range(len(A[0]) - 1):
        del M[i][j]
    return M
 
 
def det(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    signum = 1
    determinant = 0
 
    for j in range(n):
        determinant += A[0][j] * signum * det(minor(A, 0, j))
        signum *= -1
    return determinant
 
def transpose(array):
    res = []
    n = len(array)
    m = len(array[0])
    for j in range(m):
        tmp=[]
        for i in range(n):
            tmp=tmp+[array[i][j]]
        res = res+[tmp]
    return res
 
def inverse_t(A):
    n = len(A)
    result = [[0 for row in range(n)] for col in range(n)]
    for i in range(n):
        for j in range(n):
            tmp = minor(A, i, j)
            if i +j % 2 == 1:
                result[i][j] = -1 * det(tmp) / det(A)
            else:
                result[i][j] = 1 * det(tmp) / det(A)
    return transpose(result)

#=================================================================

def bin_search(a, b, epsilon, f):
    while b - a > epsilon:
        c = (a + b) / 2
        if f(b) * f(c) < 0:
            a = c
        else:
            b = c
    return (a + b) / 2

def newton_tangents(x0, epsilon, f, f_der):
    xn = x0
    xnm1 = x0-1
    q = -1
    
    for i in range(1000000):
        xnp1 = xn - f(xn) / f_der(xn)
        q = (xnp1 - xn) / (xn - xnm1)
        if math.fabs(xnp1 - xn) < epsilon:
            break
        xnm1 = xn
        xn = xnp1
    
    q = 1 / (1 - q)

    return [xn, q]

def newton_chords(x0, x1, epsilon, f):
    xn = x0
    xnm1 = x1
    q = -1
    
    for i in range(1000000):
        xnp1 = xn - f(xn) * (xn - xnm1) / (f(xn) - f(xnm1))
        q = (xnp1 - xn) / (xn - xnm1)
        if math.fabs(xnp1 - xn) < epsilon:
            break
        xnm1 = xn
        xn = xnp1
        
    q = 1 / (1 - q)
    
    return xn#[xn, q]

def newton_tangent(x0, epsilon, f, f_der):
    xn = x0
    xnm1 = x0 - 1
    der = f_der(x0)
    q = -1
    
    for i in range(1000000):
        xnp1 = xn - f(xn) / der
        q = (xnp1 - xn) / (xn - xnm1)
        if math.fabs(xnp1 - xn) < epsilon:
            break
        xnm1 = xn
        xn = xnp1
        
    q = 1 / (1 - q)
    
    return xn#[xn, q]

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

def newton_multidimensional(F, Jacob, X0, epsilon):
    Xn = veccopy(X0)
    
    for i in range(1000000):
        # JacobXn * Xk = Jacob_Xn*Xk - F_Xn
        # A*X = B
        # A = vectorsub
        # B = Jacob_Xn*Xk - F_Xn
        Jacob_Xn = matrix_from_multidimensional_function(Jacob, Xn)
        F_Xn = vector_from_multidimensional_function(F, Xn)
        
        Jacob_Xn_Xn = vectormul(Jacob_Xn, Xn)
        
        A = Jacob_Xn
        B = vectorsub(Jacob_Xn_Xn, F_Xn)
        
        gauss_simple(A, B)
        Xnp1 = B
        
        if delta(Xn, Xnp1) < epsilon:
            break
        
        Xn = Xnp1
    
    return Xn

# tests
def func(x):
    return (x-1)*(x-1)*(x-1)#math.cos(x) - x*x*x

def func_der(x):
    return 3*(x-1)*(x-1)#-math.sin(x)-3*x*x

# multi tests
def U(X):
    return (X[0]-2)*(X[0]-2) + (X[1]-4)*(X[1]-4)

def V(X):
    return 2*X[0]-X[1]

def Ux(X):
    return 2*(X[0]-2)

def Uy(X):
    return 2*(X[1]-4)

def Vx(X):
    return 2

def Vy(X):
    return -1

F = [U, V]
Jacob = [
            [Ux, Uy],
            [Vx, Vy]
        ]

X0 = [0.7, 0.18]

# Tests
print(bin_search(0, 3, 0.001, func))
print(newton_tangents(0.5, 0.001, func, func_der))
print(newton_tangent(0.5, 0.001, func, func_der))
print(newton_chords(200, 30, 0.001, func))
print(newton_multidimensional(F, Jacob, X0, 0.001))

