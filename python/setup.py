#!/usr/bin/env python

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
      packages=['getdist','getdist.gui','paramgrid'],
      scripts=['GetDist'],
      package_data={'getdist': ['getdist/analysis_defaults.ini']},
      requires=[
          'numpy',
          'matplotlib',
          "scipy (>=0.11.0)"]
 #These are optional
 #       'GUI':  ["PySide"],
 #       'fastread': ["pandas>=0.14.0"]}
      )
