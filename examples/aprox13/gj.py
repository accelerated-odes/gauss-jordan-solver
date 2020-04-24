from SparseGaussJordan import GaussJordan

GJ = GaussJordan(structure_file="maskfile", compressed_sparse_row=False, out_py="solver.py", out_f95="solver.f90", out_cpp="actual_linear_solver.H", cpp_use_arraynd=True, cse=True, smp=False)
