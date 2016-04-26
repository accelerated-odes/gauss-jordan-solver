#!/usr/bin/env python
import numpy as np
import argparse
from gauss_jordan import GaussJordan, Row

def random_sample(GJ, scale, verbose=False):
    x = scale * np.random.rand(GJ.nvars)
    A = np.zeros((GJ.nvars, GJ.nvars), dtype=np.float64)
    for j, r in enumerate(GJ.sparsity):
        for k in range(GJ.nvars):
            if not r.zero_at(k):
                frand = np.random.rand()
                while frand==0.0:
                    # Get a nonzero random number for element of A.
                    frand = np.random.rand()
                A[j][k] = scale * frand
                if verbose:
                    print('A[{}][{}] = random'.format(j, k))
            else:
                A[j][k] = 0.0
                if verbose:
                    print('A[{}][{}] = 0'.format(j, k))
    return x, A

def test_gj_solve(GJ, ntest, scale, tol):
    # Test the Gauss-Jordan solver by randomizing matrix elements and verifying output
    xresd_accumulate = np.empty((ntest, GJ.nvars), dtype=np.float64)
    xtrue, A = random_sample(GJ, scale, verbose=True) # Print which elements get set
    print(A)
    for i in range(ntest):
        xtrue, A = random_sample(GJ, scale)
        b = np.dot(A, xtrue)
        xtest = gauss_jordan_solve(A, b)
        xresd = np.absolute(xtrue-xtest) # Absolute value of difference
        xresd_accumulate[i][:] = xresd[:]
    # Get statistics
    xresd_ave = np.average(xresd_accumulate, axis=0)
    xresd_std = np.std(xresd_accumulate, axis=0)
    # Print statistics
    print('Residual Average +/- StdDev for solutions of the system Ax=b.')
    print('Using {} random iterations.'.format(ntest))
    is_okay = True
    for i, xstat in enumerate(zip(xresd_ave, xresd_std)):
        print('Residual for X[{}]: {} +/- {}'.format(i, xstat[0], xstat[1]))
        if xstat[0] > tol or xstat[1] > tol:
            is_okay = False
        if not is_okay:
            print('Tolerance ERROR for X[{}]!'.format(i))
    if is_okay:
        print('System OK!\nAverage and standard deviation of all residuals is below the tolerance {}'.format(tol))
        
if __name__=='__main__':
    ### MAIN ###
    parser = argparse.ArgumentParser()
    parser.add_argument('structurefile', type=str,
                        help="""Name of the file specifying sparsity pattern 
                        of the matrix on which to test the gauss-jordan solution.""")
    parser.add_argument('-n', type=int, default=15,
                        help='Number of random matrices to test. (Default is 15).')
    parser.add_argument('-scale', type=float, default=100.0,
                        help='Scaling factor with which to multiply the random matrix elements.')
    parser.add_argument('-tol', type=float, default=1.0e-12,
                        help="""Round-off error tolerance for residuals (both average and stdev) 
                        below which the solution is OK. (Default is 1.0E-12).""") 
    parser.add_argument('-smp', action='store_true',
                        help="""Attempt to simplify solution. Can be very slow, but if possible, 
                        will reduce the number of operations required for the solution.""")
    parser.add_argument('-cse', action='store_true',
                        help="""Execute Common Subexpression Elimination. 
                        (After simplification if the -smp option is present.) 
                        This may be more helpful for Python than a compiled language 
                        with a compiler capable of weighing the operation and memory costs.""")
    args = parser.parse_args()

    # Generate the python script file 'test.py' to solve the system
    sfile = args.structurefile
    GJ = GaussJordan(sfile, 'test.py', 'test.f95', args.smp, args.cse)
    from test import gauss_jordan_solve
    test_gj_solve(GJ, args.n, args.scale, args.tol)
