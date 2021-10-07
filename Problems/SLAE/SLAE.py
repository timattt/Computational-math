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
    return [[int(elem) for elem in input().split()] for i in range(n)]
    
def inputVector(n):
    return [int(elem) for elem in input().split()]

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
    
    for i in range(len(B)):
        C[i] = A[i] * b
        
    return C

def dot(A, B):
    res = 0
    for i in range(len(A)):
        res += A[i] * B[i]
    return res

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
    
def jacob(A, b):
    n = len(A)
    D = identity(n)
    
    for i in range(n):
        D[i][i] = 1 / A[i][i]
        
    
    B = matrixsub(identity(n), matrixmul(D, A))
    g = vectormul(D, b)
    
    FancyPrint(B, g, None)
    
    Xk = [0 for row in range(len(A))]
    
    for i in range(10000):
        Xk1 = vectorsum(vectormul(B, Xk), g)
        if delta(Xk, Xk1) < 0.01:
            break
        Xk = Xk1
        
    for i in range(len(b)):
        b[i] = Xk[i]
        
    return 1
    
def zeidel(A, b):
    n = len(b)
    C = [[0 for row in range(n)] for col in range(n)]
    
    for i in range(n):
        for j in range(n):
            if i != j:
                C[i][j] = - A[i][j] / A[i][i]
                
    d = [0 for row in range(n)]
    for i in range(n):
        d[i] = b[i] / A[i][i]
    
    
    Xk = [0 for row in range(n)]
    for it in range(10000):
        Xk1 = [0 for row in range(n)]
        for i in range(n):
            Xk1[i] = d[i];
            for j in range(0, i):
                Xk1[i] += C[i][j] * Xk1[j]
            for j in range(i, n):
                Xk1[i] += C[i][j] * Xk[j]
        if delta(Xk, Xk1) < 0.01:
            break
        Xk = Xk1
                
    for i in range(len(b)):
        b[i] = Xk[i]
    return 1

def fastDescend(A, b):
    n = len(b)
    Tn = 0
    rn = [0 for row in range(n)]
    Xn = [0 for row in range(n)]
    for i in range(100000):
        rn = vectorsub(vectormul(A, Xn), b)
        Tn = dot(rn, rn) / dot(vectormul(A, rn), rn)
        Xn1 = vectorsub(Xn, vectormulf(rn, Tn))
        if delta(Xn, Xn1) < 0.01:
            break
        Xn = Xn1
    
    for i in range(len(b)):
        b[i] = Xn[i]
        
    return 1

def leastResidual(A, b):
    n = len(b)
    Tn = 0
    rn = [0 for row in range(n)]
    Xn = [0 for row in range(n)]
    for i in range(1000000):
        rn = vectorsub(vectormul(A, Xn), b)
        Tn = dot(vectormul(A, rn), rn) / dot(vectormul(A, rn), vectormul(A, rn))
        Xn1 = vectorsub(Xn, vectormulf(rn, Tn))
        if delta(Xn, Xn1) < 0.01:
            break
        Xn = Xn1
    
    for i in range(len(b)):
        b[i] = Xn[i]
        
    return 1

print("SLAE Solver")

print("Input n...")
n = int(input())

print("Input matrix...")
A = inputMatrix(n)

print("Input vector...")
B = inputVector(n)

print("Cjoose method to solve [1 - Gauss simple, 2 - Gauss with lead element, 3 - Jacob, 4 - Zeidel], 5 - fast descend, 6 - least residual")
m = int(input())

if m == 1:
    if gauss_simple(A, B):
        print("Result:")
        print(B)
    else:
        print("Bad matrix")
if m == 2:
    if gauss_lead(A, B):
        print("Result:")
        print(B)
    else:
        print("Bad matrix")
if m == 3:
    if jacob(A, B):
        print("Result:")
        print(B)
    else:
        print("Bad matrix")
if m == 4:
    if zeidel(A, B):
        print("Result:")
        print(B)
    else:
        print("Bad matrix")
if m == 5:
    if fastDescend(A, B):
        print("Result:")
        print(B)
    else:
        print("Bad matrix")
if m == 6:
    if leastResidual(A, B):
        print("Result:")
        print(B)
    else:
        print("Bad matrix")