import numpy as np

def gauss_jordan_solve(A, b):
    x = np.empty(4, dtype=np.float64)

    x[0] = -(-A[0][1]*(A[1][2] - A[0][2]*A[1][0]/A[0][0])/(A[0][0]*(A[1][1] - A[0][1]*A[1][0]/A[0][0])) + A[0][2]/A[0][0])*(b[2] - (A[2][1] - A[0][1]*A[2][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[2][0]*b[0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (-(-A[0][1]*(A[1][2] - A[0][2]*A[1][0]/A[0][0])/(A[0][0]*(A[1][1] - A[0][1]*A[1][0]/A[0][0])) + A[0][2]/A[0][0])*(A[2][3] - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[2][0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - A[0][1]*(A[1][3] - A[0][3]*A[1][0]/A[0][0])/(A[0][0]*(A[1][1] - A[0][1]*A[1][0]/A[0][0])) + A[0][3]/A[0][0])*(b[3] - (A[3][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[3][0]/A[0][0])*(b[2] - (A[2][1] - A[0][1]*A[2][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[2][0]*b[0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (A[3][1] - A[0][1]*A[3][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[3][0]*b[0]/A[0][0])/(A[3][3] - (A[2][3] - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[2][0]/A[0][0])*(A[3][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[3][0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[3][0]/A[0][0]) - A[0][1]*(b[1] - A[1][0]*b[0]/A[0][0])/(A[0][0]*(A[1][1] - A[0][1]*A[1][0]/A[0][0])) + b[0]/A[0][0]
    x[1] = -(-(A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][3] - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[2][0]/A[0][0])/((A[1][1] - A[0][1]*A[1][0]/A[0][0])*(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0])) + (A[1][3] - A[0][3]*A[1][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]))*(b[3] - (A[3][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[3][0]/A[0][0])*(b[2] - (A[2][1] - A[0][1]*A[2][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[2][0]*b[0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (A[3][1] - A[0][1]*A[3][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[3][0]*b[0]/A[0][0])/(A[3][3] - (A[2][3] - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[2][0]/A[0][0])*(A[3][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[3][0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[3][0]/A[0][0]) - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(b[2] - (A[2][1] - A[0][1]*A[2][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[2][0]*b[0]/A[0][0])/((A[1][1] - A[0][1]*A[1][0]/A[0][0])*(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0])) + (b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0])
    x[2] = -(A[2][3] - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[2][0]/A[0][0])*(b[3] - (A[3][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[3][0]/A[0][0])*(b[2] - (A[2][1] - A[0][1]*A[2][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[2][0]*b[0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (A[3][1] - A[0][1]*A[3][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[3][0]*b[0]/A[0][0])/((A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0])*(A[3][3] - (A[2][3] - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[2][0]/A[0][0])*(A[3][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[3][0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[3][0]/A[0][0])) + (b[2] - (A[2][1] - A[0][1]*A[2][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[2][0]*b[0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0])
    x[3] = (b[3] - (A[3][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[3][0]/A[0][0])*(b[2] - (A[2][1] - A[0][1]*A[2][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[2][0]*b[0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (A[3][1] - A[0][1]*A[3][0]/A[0][0])*(b[1] - A[1][0]*b[0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[3][0]*b[0]/A[0][0])/(A[3][3] - (A[2][3] - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[2][0]/A[0][0])*(A[3][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[3][0]/A[0][0])/(A[2][2] - (A[1][2] - A[0][2]*A[1][0]/A[0][0])*(A[2][1] - A[0][1]*A[2][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][2]*A[2][0]/A[0][0]) - (A[1][3] - A[0][3]*A[1][0]/A[0][0])*(A[3][1] - A[0][1]*A[3][0]/A[0][0])/(A[1][1] - A[0][1]*A[1][0]/A[0][0]) - A[0][3]*A[3][0]/A[0][0])

    return x