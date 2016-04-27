(time ./test_gj_solve.gfortran.exe structure_4x4-8 1000000 1 1) > gauss_jordan_1.time 2>&1

(time ./test_gj_solve.gfortran.exe structure_4x4-8 1000000 1 2) > gauss_jordan_2.time 2>&1

(time ./test_gj_solve.gfortran.exe structure_4x4-8 1000000 2 1) > lapack_dgesv_1.time 2>&1

(time ./test_gj_solve.gfortran.exe structure_4x4-8 1000000 2 2) > lapack_dgesv_2.time 2>&1
