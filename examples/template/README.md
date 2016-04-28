# Gauss-Jordan Template

This template is useful for generating solvers for new matrix sparsity
masks for the matrix A in the system Ax=b.

See ```gj.py``` for a description of the options available and to
change any output or algebraic processing settings.

## USING:

* Edit ```maskfile``` to contain the array sparsity mask for matrix A.

Replace:

```
[number of rows]
[row 1 sparsity mask]
[row 2 sparsity mask]
...
```

With: (e.g. for a sample 4x4 matrix)

```
4
1 0 1 0
0 0 1 1
0 1 0 1
1 1 0 0
```

* Edit ```gj.py``` to change python or fortran output filenames,
  enable algebraic simplification, enable common subexpression
  elimination, or enable verbose output.

* Run ```$ python gj.py``` to create the solver output file(s)


## TESTING/TIMING:

* Run ```$ make``` to build the Fortran test program

* Run ```$ bash runtiming.sh``` to run a timing suite comparing the
  gauss-jordan solver with LAPACK's DGESV.

The test suite will produce four ```*.time``` files, two each for
the gauss-jordan solver and the DGESV solver. For each solver, the
test program will solve a set (by default 10^6) of random matrix
equations with the sparsity specified in ```maskfile```.

Testing has 2 ways of solving a set of random matrix equations:

* (1) Solve 10^6 different random systems (output in ```*_1.time``` files)

* (2) Solve the same random system 10^6 times (output in ```*_2.time``` files)

In each case, the average residual vector between the true solution
and the solver's solution is written to the output files as well as
(if applicable) the standard deviation across the solutions. The
timing reported by the GNU ```time``` utility is also included.