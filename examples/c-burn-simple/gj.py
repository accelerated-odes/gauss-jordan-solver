from gauss_jordan import GaussJordan

GJ = GaussJordan(structure_file = 'sparsity.dat',
                 out_f95 = 'gauss_jordan_module.f90',
                 smp = False,
                 cse = True,
                 verbose = False)
