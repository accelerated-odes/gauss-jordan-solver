import numpy as np

def gauss_jordan_solve(A, b):
    x = np.empty(2, dtype=np.float64)

    x[0] = -A[0][1]*(b[1] - A[1][0]*b[0]/A[0][0])/(A[0][0]*(A[1][1] - A[0][1]*A[1][0]/A[0][0])) + b[0]/A[0][0]
    x[1] = (b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0])

    return x
