from SparseGaussJordan import GaussJordan

GJ = GaussJordan(structure_file="maskfile", compressed_sparse_row=True,
                 out_cpp="solver.cpp", cpp_template="solver_template.cpp",
                 cse=True, smp=False)
