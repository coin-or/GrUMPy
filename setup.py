#!/usr/bin/env python

from setuptools import setup
import setuptools

setup(name='coinor.grumpy',
      version='0.8.4',
      description='Graphics for Understanding Mathematical Programming (GrUMPy)',
      long_description='''GrUMPy is a class for visualizing various algorithm used in solving discrete optimization problem. It has a class for dynamically generating and visualizing branch-and-bound trees that is derived from the GiMPy graph class. Using the branch-and-bound class, a user can visualize the branch-and-bound process in a number of different ways either by building the tree dynamically through direct calls to Python from the solver or by piping the output of an instrumented solver to GrUMPy for parsing. The branch-and-bound class also includes a pure Python implementation of branch and bound that is targeted at educational use.

In addition, GrUMPy includes a class for visualizing 2-dimensional polyhedra that can be used in combination with a pure Python implementation of the Gomory cutting plane algorithm to geometrically visualize the process of solving an integer program by a cutting plane algorithm. In future releases, the cutting plane visualization will be joined together with the branch-and-bound implementation to yield a full-blown visualization of the branch-and-cut algorithm.

Documentation for the API is here:

https://tkralphs.github.io/GrUMPy
''',
      author='Aykut Bulut, Ted Ralphs',
      author_email='{aykut,ted}@lehigh.edu',
      license='Eclipse Public License',
      url='https://github.com/tkralphs/GrUMPy/',
      namespace_packages=['coinor'],
      packages=['coinor.grumpy','coinor.grumpy.examples','coinor'],
      package_dir = {'coinor': 'src'},
      install_requires=['coinor.gimpy>=1.3.0', 'pulp']
     )
