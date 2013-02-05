import GetDistPlots
from pylab import *
import sys

g=GetDistPlots.GetDistPlotter('main/plot_data')
g.settings.setWithSubplotSize(3)
g.settings.axes_fontsize +=2
g.settings.lab_fontsize +=6
g.settings.solid_colors = ['#009966', '#000866', '#000866','#336600', '#006633' , 'g', 'm', 'r'] 

bases=['omegabh2','omegach2','ns','H0','sigma8']
params=['omegak','mnu','nnu','yhe','nrun','r']
defs=[0,0.06,3.046,0.247816,0,0]
bfs=[0.02209,0.1194,0.962,67.5,0.827]


roots=[['base_'+param+'_planck_CAMspec_lowl_lowLike','base_'+param+'_planck_CAMspec_lowl_lowLike_post_BAO','base_'+param+'_planck_CAMspec_lowl_lowLike_BAO'] for param in params]

g.rectangle_plot(bases,params,roots,ymarkers=defs, xmarkers=bfs)
g.export('outputs/grid_1paramext_planck_CAMspec_lowl_lowLike.eps')

