import planckStyle as s
from pylab import *
import sys

g=s.plotter
0
g.settings.setWithSubplotSize(2)
g.settings.axes_fontsize +=2
g.settings.lab_fontsize +=5
g.settings.solid_colors = ['#009966', '#000866', '#000866'] 

bases=['omegabh2','omegach2','ns','H0','sigma8']
params=['omegak','mnu','nnu','yhe','nrun','r','w']
defs=[0,0.06,3.046,0.247816,0,0,-1]
bfs=[0.02205,0.1199,0.9619,67.04,0.8347]


roots=[['base_'+param+'_planck_lowl_lowLike','base_'+param+'_planck_lowl_lowLike_post_BAO','base_'+param+'_planck_lowl_lowLike_BAO'] for param in params]

g.rectangle_plot(bases,params,roots,ymarkers=defs, xmarkers=bfs)
g.export('grid_1paramext_planck_lowl_lowLike')

