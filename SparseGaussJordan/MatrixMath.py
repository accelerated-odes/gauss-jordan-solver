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

class MatrixMath(object):
    @staticmethod
    def csr_matrix_vector_dot(GJ, A, x):
        b = np.zeros_like(x)

        icsr = 0
        for irow, r in enumerate(GJ.sparsity):
            for icol in range(GJ.nvars):
                if not r.zero_at(icol):
                    b[irow] += A[icsr] * x[icol]
                    icsr += 1
        return b

    @staticmethod
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
