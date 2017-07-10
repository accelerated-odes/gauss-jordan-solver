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
# Solve Ax=b with gauss-jordan elimination
# Requires Python 3.
import sympy
import re
import argparse

class Element2D(object):
    def __init__(self, symbase, i, j):
        self.symbase = symbase # Array name (e.g. A)
        self.i = i # 1-based row index into the 2D array
        self.j = j # 1-based col index into the 2D array
        self.pyrep = self.pythonify()
        self.fnrep = self.fortranify()

    def pythonify(self):
        return '{}[{}][{}]'.format(self.symbase, self.i-1, self.j-1)

    def fortranify(self):
        return '{}({},{})'.format(self.symbase, self.i, self.j)

class Element1D(object):
    def __init__(self, symbase, i):
        self.symbase = symbase # Array name (e.g. A)
        self.i = i # 1-based row index into the 1D array
        self.pyrep = self.pythonify()
        self.fnrep = self.fortranify()

    def pythonify(self):
        return '{}[{}]'.format(self.symbase, self.i-1)

    def fortranify(self):
        return '{}({})'.format(self.symbase, self.i)

class Row(object):
    def __init__(self, elist=None):
        self.elements = elist

    def __iter__(self):
        for es in self.elements:
            yield es
        
    def __add__(self, other):
        if isinstance(other, Row):
            if len(self.elements) != len(other.elements):
                print('ERROR: Unequal row lengths!')
                exit()
            else:
                sum_elements = []
                for es, eo in zip(self.elements, other.elements):
                    sum_elements.append(es+eo)
                return Row(sum_elements)
        elif isinstance(other, (int, float) + tuple(sympy.core.all_classes)):
            sum_elements = []
            for es in self.elements:
                sum_elements.append(es+other)
            return Row(sum_elements)
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Row):
            if len(self.elements) != len(other.elements):
                print('ERROR: Unequal row lengths!')
                exit()
            else:
                sum_elements = []
                for es, eo in zip(self.elements, other.elements):
                    sum_elements.append(es-eo)
                return Row(sum_elements)
        elif isinstance(other, (int, float) + tuple(sympy.core.all_classes)):
            sum_elements = []
            for es in self.elements:
                sum_elements.append(es-other)
            return Row(sum_elements)
        else:
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        if isinstance(other, Row):
            if len(self.elements) != len(other.elements):
                print('ERROR: Unequal row lengths!')
                exit()
            else:
                sum_elements = []
                for es, eo in zip(self.elements, other.elements):
                    sum_elements.append(eo-es)
                return Row(sum_elements)
        elif isinstance(other, (int, float) + tuple(sympy.core.all_classes)):
            sum_elements = []
            for es in self.elements:
                sum_elements.append(other-es)
            return Row(sum_elements)
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, Row):
            # Perform element-wise multiplication
            if len(self.elements) != len(other.elements):
                print('ERROR: Unequal row lengths!')
                exit()
            else:
                mul_elements = []
                for es, eo in zip(self.elements, other.elements):
                    mul_elements.append(eo*es)
                return Row(mul_elements)
        elif isinstance(other, (int, float) + tuple(sympy.core.all_classes)):
            mul_elements = []
            for es in self.elements:
                mul_elements.append(es*other)
            return Row(mul_elements)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Row):
            # Perform element-wise division
            if len(self.elements) != len(other.elements):
                print('ERROR: Unequal row lengths!')
                exit()
            else:
                div_elements = []
                for es, eo in zip(self.elements, other.elements):
                    div_elements.append(es/eo)
                return Row(div_elements)
        elif isinstance(other, (int, float) + tuple(sympy.core.all_classes)):
            div_elements = []
            for es in self.elements:
                div_elements.append(es/other)
            return Row(div_elements)
        else:
            return NotImplemented
        
    def __rtruediv__(self, other):
        if isinstance(other, Row):
            # Perform element-wise division
            if len(self.elements) != len(other.elements):
                print('ERROR: Unequal row lengths!')
                exit()
            else:
                div_elements = []
                for es, eo in zip(self.elements, other.elements):
                    div_elements.append(eo/es)
                return Row(div_elements)
        elif isinstance(other, (int, float) + tuple(sympy.core.all_classes)):
            div_elements = []
            for es in self.elements:
                div_elements.append(other/es)
            return Row(div_elements)
        else:
            return NotImplemented

    def fnzero(self):
        # Finds the first nonzero element in the Row: x
        # Returns the tuple (i, x) where i is the index of x.
        for n, e in enumerate(self.elements):
            if e!=0:
                return (n, e)
        return (-1, 0) # Return (-1, 0) if there are no nonzero elements in the Row

    def zero_at(self, i):
        # Returns True if element in ith position is zero, False otherwise
        if self.elements[i]==0:
            return True
        else:
            return False

class GaussJordan(object):
    def __init__(self, structure_file=None, out_py=None, out_f95=None, smp=None, cse=None, verbose=None):
        self.infile = structure_file
        self.cse_rep = None
        self.verbose = verbose
        self.sparsity = None # Original matrix system sparsity pattern
        self.symtab = {}

        if self.infile:
            self.readfile()
            if self.verbose:
                self.printEquation()
            self.pivotSolve()
            if smp:
                print('Simplifying')
                self.simplify()
            if cse:
                print('Eliminating CSE')
                self.elim_cse()
            if self.verbose:
                self.printSolution()
            if out_py:
                self.writeCode_Python(out_py)
            if out_f95:
                self.writeCode_Fortran95(out_f95)
                
    def readfile(self):
        # Read in the array mask and store in amat as either 'False'
        # or the string 'amat(i,j)' where i is the row number and
        # j the column number (starting at 1).
        try:
            file = open(self.infile,'r')
        except:
            raise

        lines = []
        for l in file:
            lines.append(l.strip())
        file.close()

        amat_t = []
        ncols = None
        for iline, line in enumerate(lines):
            if iline == 0:
                continue # This is the size of the square matrix, useful for Fortran
            if line=='':
                continue
            elif line[0]=='#':
                continue
            cols = []
            for e in line.split():
                if e=='0':
                    cols.append(False)
                else:
                    cols.append(True)
            if not ncols:
                ncols = len(cols)
            else:
                if ncols != len(cols):
                    print('Error: inconsistent number of columns!')
                    print(cols)
                    print(ncols)
                    exit()
            amat_t.append(cols)

        if self.verbose:
            print(amat_t)
            print(ncols)

        if len(amat_t) != ncols:
            print('Error: number of rows not equal to number of columns.')
            exit()

        self.asym = []
        for i, r in enumerate(amat_t):
            r_s = []
            for j, a_ij in enumerate(r):
                if self.verbose:
                    print(a_ij)
                if a_ij:
                    symrep = 'A_'+str(i+1)+'_'+str(j+1)+'_'
                    self.symtab[symrep] = Element2D('A', i+1, j+1)
                    r_s.append(sympy.symbols(symrep))
                else:
                    r_s.append(0)
            symrep = 'b_'+str(i+1)+'_'
            self.symtab[symrep] = Element1D('b', i+1)
            r_s.append(sympy.symbols(symrep))
            self.asym.append(Row(r_s))

        # Save a copy of the original sparsity
        self.sparsity = [r for r in self.asym]
        
    def printEquation(self):
        print('Symbolic array equation [A|b]:')
        for r in self.asym:
            s = ''
            for ri in r:
                s = s + '{:20}'.format(str(ri))
            print(s)

    def pivotSolve(self):
        self.nvars = len(self.asym)
        vpivot = 0
        # Check if already in reduced-row-echelon (minus the row swapping)
        # Pivot for variable vpivot
        for i in range(vpivot, self.nvars):
            vrows = []
            for j, r in enumerate(self.asym):
                # Make a list of rows with nonzero elements in the ith position
                if self.verbose:
                    print('Checking row {}'.format(j))
                    print('row {}, index {}: {}'.format(j, i, str(r.elements[i])))
                if not r.zero_at(i):
                    vrows.append(j)
            for j, r in enumerate(self.asym):
                n, e = r.fnzero()
                if n==i:
                    ## Pivot using jth Row
                    # Copy the matrix
                    asym_pivot = [r for r in self.asym]
                    # Normalize element in the pivot position
                    asym_pivot[j] = asym_pivot[j]/e
                    # Pivot other rows with a nonzero element in ith position
                    for vj in vrows:
                        if vj != j:
                            # Get elimination prefactor
                            prefactor = asym_pivot[vj].elements[i]
                            # Eliminate into this row
                            asym_pivot[vj] = asym_pivot[vj] - prefactor * asym_pivot[j]
                    ## The following uses only the first available elimination
                    # Replace the original matrix with the pivot matrix
                    self.asym = [r for r in asym_pivot]
                    break

        # Rearrange rows into RREF
        asym_temp = []
        for i in range(len(self.asym)):
            for r in self.asym:
                n, e = r.fnzero()
                if n == i:
                    asym_temp.append(r)
        self.asym = [r for r in asym_temp]

        # Find solution vector
        self.solution = [r.elements[-1] for r in self.asym]

    def simplify(self):
        # Simplify solutions using simplify
        # ratio=1 requires the number of operations in the simplified
        # expression to be at most the number of operations in the original
        if self.verbose:
            print('Finding simplifications...')
        for i, s in enumerate(self.solution):
            sres = sympy.simplify(s, ratio=1)
            self.solution[i] = sres
            if self.verbose:
                print('----------------------------------------')
                print('Simplify Solution {} Report:'.format(i))
                if str(sres) != str(s):
                    print('Shorter Representation Found')
                    print('Un-Simplified:')
                    print(s)
                    print('Simplified:')
                    print(sres)
                else:
                    print('No Shorter Representation Found')
                print('----------------------------------------')
                    

    def elim_cse(self):
        ## Apply common sub-expressions elimination to solution
        scratch_sym = sympy.utilities.numbered_symbols('scratch_')
        cse_rep, cse_sol = sympy.cse(self.solution, symbols=scratch_sym, order='none')
        self.cse_rep = cse_rep
        if self.verbose:
            print('--> Eliminated Common Sub-Expressions')
            print('Scratch Variables:')
            for rep in cse_rep:
                print('{} = {}'.format(rep[0], rep[1]))
            print('Solution:')
        self.solution = cse_sol

    def printSolution(self):
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('FINAL SOLUTION')
        print('')
        for i, sol in enumerate(self.solution):
            print('X{} = {}'.format(i+1, sol))
            print('')

    def pythonify(self,s):
        # Python-ify string s by replacing A_i_j_ with A[i][j]
        # and b_i_ with b[i].
        # Also start indexing at 0 instead of 1.
        sout = str(s)
        for k in self.symtab.keys():
            sout = sout.replace(k, self.symtab[k].pyrep)
        return sout

    def fortranify(self,s):
        # Fortran-ify string s by replacing A_i_j_ with A(i,j)
        # and b_i_ with b(i)
        sout = str(s)
        for k in self.symtab.keys():
            sout = sout.replace(k, self.symtab[k].fnrep)
        return sout
    
    def writeCode_Python(self, outname):
        ## Write Python Solver Function
        try:
            fo = open(outname, 'w')
        except:
            raise
        indent = ' '*4
        fo.write('import numpy as np\n')
        fo.write('\n')
        fo.write('def gauss_jordan_solve(A, b):\n')
        fo.write('{}x = np.empty({}, dtype=np.float64)\n'.format(indent, self.nvars))
        fo.write('\n')
        if self.cse_rep:
            for rep in self.cse_rep:
                fo.write('{}{} = {}\n'.format(indent, rep[0], self.pythonify(rep[1])))
            fo.write('\n')
        for i, sol in enumerate(self.solution):
            fo.write('{}x[{}] = {}\n'.format(indent, i, self.pythonify(sol)))
        fo.write('\n')
        fo.write('{}return x\n'.format(indent))
        fo.close()

    def writeCode_Fortran95(self, outname):
        ## Write Fortran-95 Solver Subroutine
        try:
            fo = open(outname, 'w')
        except:
            raise
        indent = ' '*2
        fo.write('module gauss_jordan_module\n')
        fo.write('{}implicit none\n'.format(indent))
        fo.write('\n')
        fo.write('contains\n')
        fo.write('\n')
        fo.write('{}subroutine gauss_jordan_solve(A, x, b)\n'.format(indent))
        fo.write('{}double precision, dimension({},{}), intent(in) :: A\n'.format(
            indent*2, self.nvars, self.nvars))
        fo.write('{}double precision, dimension({}), intent(out) :: x\n'.format(
            indent*2, self.nvars))
        fo.write('{}double precision, dimension({}), intent(in) :: b\n'.format(
            indent*2, self.nvars))
        if self.cse_rep:
            for rep in self.cse_rep:
                fo.write('{}double precision :: {}\n'.format(indent*2, rep[0]))
            fo.write('\n')
            for rep in self.cse_rep:
                rep_value = self.fortranify(sympy.fcode(rep[1], precision = 15,
                                                        source_format = 'free',
                                                        standard = 95))
                fo.write('{}{} = {}\n'.format(indent*2, rep[0], rep_value))
        fo.write('\n')
        for i, sol in enumerate(self.solution):
            sol_value = self.fortranify(sympy.fcode(sol, precision = 15,
                                                    source_format = 'free',
                                                    standard = 95))
            fo.write('{}x({}) = {}\n'.format(indent*2, i+1, sol_value))
        fo.write('{}end subroutine gauss_jordan_solve\n'.format(indent))
        fo.write('end module gauss_jordan_module\n')
        fo.close()

## Run
if __name__=='__main__':
    ## Read CMD Line Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('maskfile', type=str,
                        help='Name of the input mask file for the array A in A*x=b.')
    parser.add_argument('-py', type=str,
                        help='Name of the Python output file to generate.')
    parser.add_argument('-f95', type=str,
                        help='Name of the Fortran-95 output file to generate.')
    parser.add_argument('-smp', action='store_true',
                        help='Attempt to simplify solution. Can be very slow, but if possible, will reduce the number of operations required for the solution.')
    parser.add_argument('-cse', action='store_true',
                        help='Execute Common Subexpression Elimination. (After simplification if the -smp option is present.) This is pretty fast.')
    parser.add_argument('-v', action='store_true',
                        help='Enable verbose output.')
    args = parser.parse_args()
    
    if args.maskfile:
        GJ = GaussJordan(args.maskfile, args.py, args.f95, args.smp, args.cse, args.v)
