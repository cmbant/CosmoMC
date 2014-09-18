import pandas as pd
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

rootname = 'WL_result_2D_v4.dat'
tolerance = 2.0

cols = ['theta1','theta2','lin','nl1','nl2','nl3','nl4']
data = pd.read_csv(rootname,names=cols)
lin = data[(data.theta1<200.0)].pivot(index='theta2',columns='theta1',values='lin')
nl1 = data[(data.theta1<200.0)].pivot(index='theta2',columns='theta1',values='nl1')
nl2 = data[(data.theta1<200.0)].pivot(index='theta2',columns='theta1',values='nl2')
nl3 = data[(data.theta1<200.0)].pivot(index='theta2',columns='theta1',values='nl3')
nl4 = data[(data.theta1<200.0)].pivot(index='theta2',columns='theta1',values='nl4')

baseline = nl4
plots = [lin,nl1,nl2,nl3]
titles = ['HL4 - Lin','HL4 - HL1','HL4 - HL2','HL4 - HL3']

fig = plt.figure()

plot_log = False

for i in range(0,len(plots)):
    ax = fig.add_subplot(2,2,i+1)
    x = baseline.columns
    y = baseline.index
    if plot_log:
        z = np.log10(np.abs(baseline.values-plots[i].values))
        v = np.linspace(0, 2, 50, endpoint=True)
        cf = plt.contourf(x,y,z,v,cmap=plt.cm.Oranges)
    else:
        z = baseline.values-plots[i].values
        v = np.linspace(-tolerance, tolerance, 100, endpoint=True)
        cf = plt.contourf(x,y,z,v,cmap=plt.cm.jet)
    cbar = plt.colorbar(cf)
    ax.set_xlabel(r'$\theta^{\plus}_c$')
    ax.set_ylabel(r'$\theta^{\minus}_c$')
    ax.set_xlim([0,100])
    ax.set_ylim([0,200])
    ax.tick_params(axis='both',which='major',labelsize=10)
    ax.set_title(titles[i],fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    plt.plot(np.array([17.0,17.0]),np.array([0.0,1000.0]),linestyle='--',color='k')
    plt.plot(np.array([0.0,1000.0]),np.array([53.0,53.0]),linestyle='--',color='k')

plt.show()
