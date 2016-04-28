#!/usr/bin/env python
"""
Copyright (c) 2015, Michael Zingale 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
# automatically generate Makefile dependencies for Fortran 90 source.
#
# this will output all the dependency pairs amongst the source files.
#
# M. Zingale (2012-03-21)

import sys
import re
import string
import os
import argparse

def doit(prefix, program, files):

    # regular expression for ' use modulename, only: stuff, other stuff'
    # see (txt2re.com)
    use_re = re.compile("( *)(use)(\s+)((?:[a-z_][a-z_0-9]+))", 
                        re.IGNORECASE|re.DOTALL)

    module_re = re.compile("( *)(module)(\s+)((?:[a-z][a-z_0-9]+))",
                           re.IGNORECASE|re.DOTALL)

    module_proc_re = re.compile("( *)(module)(\s+)(procedure)(\s+)((?:[a-z][a-z_0-9]+))",
                                re.IGNORECASE|re.DOTALL)
    
    program_re = re.compile("( *)(program)(\s+)((?:[a-z_][a-z_0-9]+))", 
                        re.IGNORECASE|re.DOTALL)

    # first parse the files and find all the module statements.  Keep a
    # dictionary of 'module name':filename.
    modulefiles = {}
    excludefiles = []
    for file in files:

        f = open(file, "r")
        
        for line in f:

            # strip off the comments
            idx = line.find("!")
            line = line[:idx]

            rebreak = module_re.search(line)
            rebreak2 = module_proc_re.search(line)
            if rebreak and not rebreak2:
                modulefiles[rebreak.group(4)] = file
            reprogram = program_re.search(line)
            if reprogram and reprogram.group(4)!=program:
                excludefiles.append(file)
        f.close()

    # go back through the files now and look for the use statements.
    # Assume only one use statement per line.  Ignore any only clauses.
    # Build a list of dependencies for the current file and output it.
    f90sources = []
    for file in files:
        if file in excludefiles:
            continue

        f90sources.append(os.path.basename(file))
        f = open(file, "r")

        for line in f:

            # strip off the comments
            idx = line.find("!")
            line = line[:idx]

            rebreak = use_re.search(line)
            if rebreak:
                if rebreak.group(4) not in modulefiles.keys():
                    print(prefix+os.path.basename(file).replace(".f90", ".o") + ':')
                else:
                    print(prefix+os.path.basename(file).replace(".f90", ".o") + ':' +
                          prefix+os.path.basename(modulefiles[rebreak.group(4)]).replace(".f90", ".o"))
        f.close()
        print(" ")

    print("F90SOURCES := {}".format(" ".join(f90sources)))
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        help="prefix to prepend to each dependency pair, e.g., for a build directory",
                        default="")
    parser.add_argument("--program",
                        help="name of program for which to prepare compilation dependencies",
                        default="")
    parser.add_argument("files", metavar="source files", type=str, nargs="*",
                        help="F90 source files to find dependencies amongst")

    args = parser.parse_args()
    doit(args.prefix, args.program, args.files)



