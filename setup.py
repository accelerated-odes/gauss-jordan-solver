from setuptools import setup, find_packages

setup(name='gauss-jordan-solver',
      version='1.2.0',
      description='Python code generator for a sparse gauss jordan linear solver',
      url='https://github.com/dwillcox/gauss-jordan-solver',
      author='Donald Willcox',
      author_email='eugene.willcox@gmail.com',
      license='BSD',
      packages=find_packages(),
      package_data={"gauss-jordan-solver": ["SparseGaussJordan/*", "util/*"]},
      install_requires=['numpy', 'sympy'],
      zip_safe=False)
