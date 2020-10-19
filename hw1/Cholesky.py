def solve_equation(A, b):
    '''
    return the solution x of the equation Ax=b (A is required to be positive definite, Hermitian)
    '''
    H_d = cholesky_decomposition(A)

    n = len(b)
    y = [0] * n
    for j in range(0, n):
        y[j] = (b[j] - sum([y[k] * H_d[j][k] for k in range(0, j)])) / H_d[j][j]
    #print(y)

    x = [0] * n
    for j in range(n-1, -1, -1):
        x[j] = (y[j] - sum([x[k] * H_d[k][j] for k in range(n-1, j, -1)])) / H_d[j][j] # change the two indexs of H\dagger
    return x


def cholesky_decomposition(A):
    '''
    do Cholesky decomposition to A, return H\dagger as H_d
    '''
    n = len(A)
    H_d = [[0] * n for p in range(n)]
    H_d[0][0] = A[0][0] ** 0.5
    for i in range(1, n):
        for j in range(0, i):
            H_d[i][j] = (A[i][j] - sum([H_d[i][k] * H_d[j][k] for k in range(0, j)])) / H_d[j][j]
        H_d[i][i] = (A[i][i] - sum([H_d[i][k] ** 2 for k in range(0, i)])) ** 0.5
    return H_d

if __name__ == '__main__':
    print(cholesky_decomposition([[1, 1/2, 1/3], [1/2, 1/3, 1/4], [1/3, 1/4, 1/5]]))
