from SparseGaussJordan import GaussJordan

GJ = GaussJordan(structure_file="maskfile", compressed_sparse_row=False, out_py="solver.py", out_f95="solver.f90", cse=True, smp=False)
