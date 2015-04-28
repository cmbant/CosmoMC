#!/usr/bin/env python

from distutils.core import setup

setup(name='GetDist',
      version='1.1',
      description='GetDist Monte Carlo sample analysis, plotting and GUI',
      author='Antony Lewis',
      url="https://cosmologist.info/cosmomc/",
      packages=['getdist','paramgrid'],
      scripts=['GetDist'],
      install_requires=[
          'numpy',
          'matplotlib',
          "scipy >= 0.11.0",
      ],
      )
