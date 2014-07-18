import planckStyle as s
from pylab import *

g=s.getSinglePlotter()

g.make_figure(1, xstretch=1.3)
roots=['base_nnu_meffsterile_'+s.defdata_TT]

s8 =  np.arange(0.5, 1,0.01)


def PLSZ(sigma):
    #from Anna 18/7/2014
    return  ((0.757+0.013*sigma)/s8)**(1/0.3)*0.32

def CFTHlens(sigma):
    return  ((0.79+0.03*sigma)/s8)**(1/0.6)*0.27

def galaxygalaxy(sigma): #mandelbaum
    return  ((0.8+0.05*sigma)/s8)**(1/0.25)*0.25

def plotBounds(data, c):
    fill_between(s8, data(-2), data(2), facecolor=c, alpha=0.15, edgecolor=c, lw=0)
    fill_between(s8, data(-1), data(1), facecolor=c, alpha=0.25, edgecolor=c, lw=0)
    
    
#plotBounds(galaxygalaxy,'green')    
plotBounds(CFTHlens,'burlywood')    
plotBounds(PLSZ,'gray')    



g.plot_3d(roots, ['sigma8', 'omegam','H0'])

g.export()
