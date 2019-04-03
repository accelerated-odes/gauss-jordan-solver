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

* Neither the name of the copyright holder nor the names of its
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
from .SparseGaussJordan import GaussJordan
from .Row import Row
from .MatrixMath import MatrixMath

class RandomTesting(object):
    @staticmethod
    def test_gj_solve(GJ, compressed_sparse_row, ntest, scale, tol, gauss_jordan_solve, verbose_testing=False):
        # Test the Gauss-Jordan solver by randomizing matrix elements and verifying output
        xresd_accumulate = np.empty((ntest, GJ.nvars), dtype=np.float64)
        xtrue, A = MatrixMath.random_sample(GJ, compressed_sparse_row, scale, verbose=verbose_testing)

        for i in range(ntest):
            xtrue, A = MatrixMath.random_sample(GJ, compressed_sparse_row, scale)
            if compressed_sparse_row:
                b = MatrixMath.csr_matrix_vector_dot(GJ, A, xtrue)
            else:
                b = np.dot(A, xtrue)
            xtest = gauss_jordan_solve(A, b)
            xresd = np.absolute(xtrue-xtest)
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
            if xstat[0] > tol or max(abs(xstat[0]+xstat[1]),abs(xstat[0]-xstat[1])) > tol:
                is_okay = False
                print('Tolerance ERROR for X[{}]!'.format(i))
                print('Tolerance is {}'.format(tol))
        if is_okay:
            print('System OK!\nAverage and standard deviation of all residuals is below the absolute tolerance {}'.format(tol))


