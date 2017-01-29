''' Setup file for ising'''

from setuptools import setup

setup(name='ising',
      version='0.4.2',
      description='Basic 2d Ising model',
      url='https://github.com/DKarandikar/ising-model',
      author='Daniel Karandikar',
      license='MIT',
      packages=['ising'],
      install_requires=[
          'numpy',
          'matplotlib',
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      include_package_data=True,
      zip_safe=False)
