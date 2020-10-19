# GEM without fulcrum selection, which is appropriate to our problem.

def solve_equation(A, b):
    '''
    return the solution x of the equation Ax=b
    return -1 if 0 is found in diagonal
    '''
    if gaussian_elimination(A, b) == -1:
        return -1
    else:
        return back_substitution(A, b)

def gaussian_elimination(A, b):
    '''
    transform A into an upper triangular matrix
    return -1 if 0 is found in diagonal
    '''
    for k in range(0, len(b)-1):
        if A[k][k] == 0:
            print("find 0 in diagonal")
            return -1
        for i in range(k+1, len(b)):
            l = -A[i][k] / A[k][k]
            for j in range(0, len(b)):
                A[i][j] = A[i][j] + l * A[k][j]
            b[i] = b[i] + l * b[k]

def back_substitution(A, b):
    x = [0] * len(b)
    for i in range(len(b)-1, -1, -1):
        if A[i][i] == 0:
            print("find 0 in diagonal")
            return -1
        x[i] = b[i] - sum([A[i][k] * x[k] for k in range(len(b)-1, i, -1)])
        x[i] = x[i] / A[i][i]
    return x

if __name__=='__main__':
    print(solve_equation([[1, 1/2], [1/2, 1/3]], [1, 1]))
