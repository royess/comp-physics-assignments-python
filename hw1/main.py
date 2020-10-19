# using Python3
# Python 3.6.3 (v3.6.3:2c5fed86e0, Oct  3 2017, 00:32:08)
# [GCC 4.2.1 (Apple Inc. build 5666) (dot 3)] on darwin

import math
import GEM
import Cholesky

def hilbert_matrix(n):
    # generate a Hilbert matrix of nxn
    hm = [ [0] * n for p in range(n)]
    for i in range(n):
        for j in range(n):
            hm[i][j] = 1 / (i + j + 1) #indexs start from 0
    return hm

def distance(vec1, vec2):
    return (sum([(x1 - x2) ** 2 for x1, x2 in zip(vec1, vec2)])) ** 0.5

def error(solution):
    # compute the error of a solution
    n = len(solution)
    hm = hilbert_matrix(n)
    result = [sum([hm[i][k] * solution[k] for k in range(n)]) for i in range(n)] # do matrix multiplication
    return abs(distance(result, [1] * n) / n ** 0.5)

if __name__ == '__main__':
    print('Calculate the solution given by GEM and Cholesky decomposition,\
     and their relative deviations, errors.', '\n')
    for n in range(1, 21):
        # n from 1 to 20, solve the equation respectively in GEM and Cholesky decomposition
        GEM_solution = GEM.solve_equation(A = hilbert_matrix(n), b = [1]*n)
        Cholesky_solution = Cholesky.solve_equation(A = hilbert_matrix(n), b = [1]*n)

        print("n=", n)
        print("GEM:", GEM_solution, "\nerror:", error(GEM_solution), '\n')
        print("Cholesky:", Cholesky_solution,\
         "\nerror:", error(Cholesky_solution), '\n')
        print("Relative Deviation:", distance(GEM_solution,\
         Cholesky_solution) / distance([0]*n, Cholesky_solution))
        if error(GEM_solution) != 0:
            print("error of Chelosky/error of GEM:",\
             error(Cholesky_solution) / error(GEM_solution))
        print('--------------------------------')

    '''
    The code below is used to compute and storage the error of Chelosky/error of GEM
    For this computation, n varies from 3 to 200, so time cost is greater
    '''

    '''
    with open("error_ratio", mode='w+') as file: #if run in win, should add '.txt'
        maximum = 0
        minimum = 0
        postive_number = 0
        for n in range(3, 201):
            res = math.log(
            error(Cholesky.solve_equation(hilbert_matrix(n), [1]*n)) / error(GEM.solve_equation(hilbert_matrix(n), [1]*n)), 10)
            if res > 0:
                postive_number = postive_number + 1
            maximum = max(res, maximum)
            minimum = min(res, minimum)
            file.write('('+ str(n) + ','+ str(res) + ') ')
        print(postive_number, maximum, minimum)
    '''
