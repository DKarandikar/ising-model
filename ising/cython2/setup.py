''' Setup file '''

from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    name='Movement scripts',
    ext_modules=cythonize("move_cy.pyx"),
    include_dirs=[numpy.get_include()]
)