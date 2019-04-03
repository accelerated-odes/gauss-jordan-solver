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

    def cify(self):
        return self.pythonify()

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

    def cify(self):
        return self.pythonify()

