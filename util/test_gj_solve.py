#!/usr/bin/env python
"""
Copyright (c) 2016, Donald E. Willcox
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of gauss-jordan-solver nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import numpy as np
import os
import argparse
from SparseGaussJordan import GaussJordan, Row

def csr_matrix_vector_dot(GJ, A, x):
    b = np.zeros_like(x)
    
    icsr = 0
    for irow, r in enumerate(GJ.sparsity):
        for icol in range(GJ.nvars):
            if not r.zero_at(icol):
                b[irow] += A[icsr] * x[icol]
                icsr += 1
    return b

def random_sample(GJ, compressed_sparse_row, scale, verbose=False):
    x = scale * np.random.rand(GJ.nvars)
    if compressed_sparse_row:
        A = np.zeros(GJ.get_number_nonzero(), dtype=np.float64)
    else:
        A = np.zeros((GJ.nvars, GJ.nvars), dtype=np.float64)
    icsr = 0
    for j, r in enumerate(GJ.sparsity):
        for k in range(GJ.nvars):
            if not r.zero_at(k):
                frand = np.random.rand()
                while frand==0.0:
                    # Get a nonzero random number for element of A.
                    frand = np.random.rand()
                if compressed_sparse_row:
                    A[icsr] = scale * frand
                else:
                    A[j][k] = scale * frand
                if verbose:
                    print('A[{}][{}] = random'.format(j, k))
                icsr += 1
            else:
                if not compressed_sparse_row:
                    A[j][k] = 0.0
                if verbose:
                    print('A[{}][{}] = 0'.format(j, k))
    return x, A

def test_gj_solve(GJ, compressed_sparse_row, ntest, scale, tol):
    # Test the Gauss-Jordan solver by randomizing matrix elements and verifying output
    xresd_accumulate = np.empty((ntest, GJ.nvars), dtype=np.float64)
    xtrue, A = random_sample(GJ, compressed_sparse_row, scale, verbose=True) # Print which elements get set
    print(A)
    for i in range(ntest):
        xtrue, A = random_sample(GJ, compressed_sparse_row, scale)
        if compressed_sparse_row:
            b = csr_matrix_vector_dot(GJ, A, xtrue)
        else:
            b = np.dot(A, xtrue)
        xtest = gauss_jordan_solve(A, b)
        xresd = xtrue-xtest # Absolute value of difference
        xresd_accumulate[i][:] = xresd[:]
    # Get statistics
    xresd_ave = np.average(xresd_accumulate, axis=0)
    xresd_std = np.std(xresd_accumulate, axis=0)
    # Print statistics
    print('----------------------------------------')
    print('Residual Average +/- StdDev for solutions of the system Ax=b.')
    print('Using {} random iterations.'.format(ntest))
    is_okay = True
    for i, xstat in enumerate(zip(xresd_ave, xresd_std)):
        print('Residual for X[{}]: {} +/- {}'.format(i, xstat[0], xstat[1]))
        if xstat[0] > tol or xstat[1] > tol:
            is_okay = False
            print('Tolerance ERROR for X[{}]!'.format(i))
            print('Tolerance is {}'.format(tol))
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
    parser.add_argument('-csr', action='store_true',
                        help='Use compressed sparse row (CSR) matrix format.')
    parser.add_argument('-smp', action='store_true',
                        help="""Attempt to simplify solution. Can be very slow, but if possible, 
                        will reduce the number of operations required for the solution.""")
    parser.add_argument('-expand', action='store_true',
                        help="""Simplify solution by expanding it and performing cancellations.""")    
    parser.add_argument('-cse', action='store_true',
                        help="""Execute Common Subexpression Elimination. 
                        (After simplification if the -smp option is present.) 
                        This may be more helpful for Python than a compiled language 
                        with a compiler capable of weighing the operation and memory costs.""")
    parser.add_argument('-keep', action='store_true',
                        help="""Do not delete 'test.py' output solver once the test is complete.""")
    args = parser.parse_args()

    # Generate the python script file 'test.py' to solve the system
    sfile = args.structurefile
    GJ = GaussJordan(sfile, args.csr, 'test.py', None, args.smp, args.expand, args.cse)
    from test import gauss_jordan_solve
    test_gj_solve(GJ, args.csr, args.n, args.scale, args.tol)

    # Delete the 'test.py' file if user didn't decide to keep it
    if not args.keep:
        here = os.getcwd()
        os.remove(os.path.join(here, 'test.py'))
