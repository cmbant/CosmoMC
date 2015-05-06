from matplotlib import rcParams, rc

# common setup for matplotlib
params = {'backend': 'pdf',
          'axes.labelsize': 9,
          'font.size': 8,
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.major.pad': 4,
          'xtick.major.pad': 4,
          'ytick.labelsize': 8,
          'text.usetex': True,
          'font.family': 'sans-serif',
          'linewidth': 0.3,
          # free font similar to Helvetica
          'font.sans-serif': 'FreeSans'}

sfmath = r'C://Work/Planck/tauNL/python/TauNL/sfmath'
# use of Sans Serif also in math mode
rc('text.latex', preamble=r'\usepackage{' + sfmath + '}')

rcParams.update(params)

