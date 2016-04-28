# gauss-jordan-solver

Writes Fortran or Python routines for solving the linear system Ax=b
ignoring zero elements of A, given its sparsity pattern.

See the ```examples``` directory for sample outputs given different
matrix sparsities.

See the ```examples/template``` directory for a starting setup for
generating and testing a new solver.


# Usage:

The only required argument is the maskfile, all others are optional.

```
$ python gauss_jordan.py [maskfile] -py [python output] -f95 [fortran output] -smp -cse -v
```

## -py: Generating Python

The ```-py``` option will generate a python file named whatever you
pass in ```[python output]``` which contains a function
(```gauss_jordan_solve```) returning the solution vector to the system
Ax=b.

## -f95: Generating Fortran

The ```-f95``` option will generate a Fortran-95 compatible file named
whatever you pass in ```[fortran output]``` which contains a
subroutine (```gauss_jordan_solve```) which finds the solution vector
to the system Ax=b.

## -smp

If present, the ```-smp``` argument uses Sympy to attempt to
algebraically simplify the solution under the constraint that
simplification must reduce the total number of operations. Sympy can
be extremely slow for dense 8x8 matrices or larger.

## -cse

If present, the ```-cse``` argument uses Sympy to find and eliminate
common sub-expressions in the solution, introducing new "scratch"
variables so that no calculation is performed more than once. In some
cases, can reduce runtime by nearly an order of magnitude (Fortran, on
a CPU). Unlike algebraic simplification, eliminating common
subexpressions with Sympy is fast.

## -v

If present, the ```-v``` argument writes verbose output to the
terminal describing the intermediate solution and other internal
steps.


# Maskfile

The maskfile specifies the sparsity layout of the matrix A composing
the left-hand-side of the system. It should read as follows:

```
[number of rows]
[row 1 sparsity mask]
[row 2 sparsity mask]
[row 3 sparsity mask]
...
```

The sparsity is a series of 1 and 0, where 0 indicates the
corresponding element of A is identically 0 and 1 indicates it is
nonzero.

## Sample Maskfile: 4x4 System

(File ```structure_4x4```)

```
4
1 0 1 0
0 0 1 1
0 1 0 1
1 1 0 0
```

## Using Sample Maskfile: 4x4 System

Writes Python and Fortran solver subroutines using common subexpression elimination.

```
$ python gauss_jordan.py structure_4x4 -py solve_4x4.py -f95 solve_4x4.f90 -cse
```


# Testing Output Solvers

## Python

The script ```test_gj_solve.py``` is provided to test the python
solver produced by ```gauss_jordan.py```. The test script takes a
sparsity file, generates a python solver, then solves a set of random
matrix equations with the provided sparsity. It compares the true
solution with the solver solution and prints residual statistics.

For example, to test the solver for the sparsity mask in
```examples/4x4-8/structure_4x4-8```, using common subexpression
elimination and 100 random matrix systems, run the following:

```$ python test_gj_solve.py examples/4x4-8/structure_4x4-8 -n 100 -cse```

## Fortran

The program ```test_gj_solve.f90``` is provided to test the fortran
solver produced by ```gauss_jordan.py```. Instructions for its use are
described in ```examples/template/README.md```.


# Requirements

* Python 3. (Developed with Python 3.5.1).

* Sympy

* numpy (for optional testing of the python solver output)
