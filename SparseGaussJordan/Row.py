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
import sympy

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
        else:
            for i, es in enumerate(self.elements):
                self.elements[i] = es + other
            return self

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
        else:
            for i, es in enumerate(self.elements):
                self.elements[i] = es - other
            return self

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
        else:
            for i, es in enumerate(self.elements):
                self.elements[i] = other - es
            return self

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
        else:
            for i, es in enumerate(self.elements):
                self.elements[i] = es * other
            return self

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
        else:
            for i, es in enumerate(self.elements):
                self.elements[i] = es / other
            return self
        
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
        else:
            for i, es in enumerate(self.elements):
                self.elements[i] = other / es
            return self

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

    def get_number_nonzero(self):
        nnz = 0
        for e in self.elements:
            if e!=0:
                nnz += 1
        return nnz
