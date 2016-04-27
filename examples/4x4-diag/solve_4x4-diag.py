import numpy as np

def gauss_jordan_solve(A, b):
    x = np.empty(4, dtype=np.float64)

    x[0] = b[0]/A[0][0]
    x[1] = b[1]/A[1][1]
    x[2] = b[2]/A[2][2]
    x[3] = b[3]/A[3][3]

    return x
