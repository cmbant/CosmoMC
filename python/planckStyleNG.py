import os
from matplotlib import rcParams, rc

# common setup for matplotlib
params = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'font.family':'sans-serif',
          'linewidth':0.3,
          # free font similar to Helvetica
          'font.sans-serif':'FreeSans'}

sfmath = r'C://Work/Planck/tauNL/python/TauNL/sfmath'
# use of Sans Serif also in math mode
rc('text.latex', preamble=r'\usepackage{' + sfmath + '}')

rcParams.update(params)

