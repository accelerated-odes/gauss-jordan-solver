from gauss_jordan import GaussJordan

GJ = GaussJordan('maskfile',          # Name of the sparsity maskfile
                 'solve.py',          # Name of output python solver (None if unwanted)
                 'solve.f90',         # Name of output fortran solver (None if unwanted)
                 False,               # If True, attempt to simplify solution algebraically
                 False,               # If True, perform Common Subexpression Elimination
                 False)               # If True, enable verbose terminal output
