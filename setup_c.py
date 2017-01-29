''' Setup file '''

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np                           # <---- New line

ext_modules = [Extension("move_cy", [r"ising\move_cy.pyx"])]

setup(
  name='Movement scripts',
  cmdclass = {'build_ext': build_ext},
  include_dirs = [np.get_include()],         # <---- New line
  ext_modules = ext_modules
)