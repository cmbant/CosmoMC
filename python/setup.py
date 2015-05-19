#!/usr/bin/env python

from __future__ import absolute_import
from distutils.core import setup
import io
import re


def find_version():
    version_file = io.open('getdist/__init__.py').read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(name='GetDist',
      version=find_version(),
      description='GetDist Monte Carlo sample analysis, plotting and GUI',
      author='Antony Lewis',
      url="https://cosmologist.info/cosmomc/",
      packages=['getdist', 'getdist.gui', 'paramgrid'],
      scripts=['GetDist', 'GetDistGUI'],
      package_data={'getdist': ['getdist/analysis_defaults.ini']},
      requires=[
          'numpy',
          'matplotlib',
          'six',
          "scipy (>=0.11.0)",
          'PySide'],
      classifiers=[
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
      ],
      # These are optional
      # 'GUI':  ["PySide"],
      # 'fastread': ["pandas>=0.14.0"]}
      )
