import planckStyle as s
from matplotlib.pyplot import *
import getdist

g = s.getSinglePlotter(chain_dir=[getdist.default_grid_root, r'C:\Tmp\Planck\KiDs', r'C:\Tmp\Planck\2017\fsig8'])

ranges = [0.2, 0.35, 0.6, 1.05]
#ranges = [0.2, 0.55, 0.4, 1.05]

pair = ['omegam', 'sigma8']

omm = np.arange(0.05, 0.7, 0.01)
s.plotBounds(omm, s.planck_lensing)

#g.plot_2d('kids450fiducial', param_pair=pair, filled=True, lims=ranges)

samples = g.sampleAnalyser.samplesForRoot('kids450fiducial', settings={'ignore_rows':0.3, 'max_scatter_points':4000})
p=samples.getParams()
samples.filter((p.A < 2.5)*(p.A > 1.7) )

s8samples = g.sampleAnalyser.samplesForRoot('fsigma-vel-theta', settings={'ignore_rows':0.3})


if False: #testing putting in Jacobian
    # print s8samples.PCA(['omegam', 'H0', 'theta'], 'LLL', 'theta')
    def jacobian(H0, omegam, sigma8):
        #assume theta \propto omm^a*H0^b
        # sigma8^2 \propto As Omegam^(1.5) H0^(3.5)
        a = 0.15
        b = 0.4
        map =np.zeros((3,3))
        map[0,0] = omegam**a*b*H0**(b-1)    #dtheta/dH0
        map[1,0] = 2*omegam*H0  #omm/dH0
        map[2,0] = -3.5/H0   #d logA/dH0
        map[0,1] = a*omegam**(a-1)*H0**b  #d theta/d omegam
        map[1,1] =  H0**2  #d omm/d omegam
        map[1,2] =  -1.5/omegam #d logA /d omegam
        map[2,0] = 0 #dtheta/dsigma8
        map[2,1] = 0  #d omm/sigma8
        map[2,2] = 2/sigma8  #d log A /d sigma8
        return np.linalg.det(map)

    import copy
    jacsamples = copy.deepcopy(s8samples)
    p=jacsamples.getParams()
    for i, (H0, omegam, sigma8) in enumerate(zip(p.H0, p.omegam, p.sigma8)):
        jacsamples.weights[i] *= jacobian(H0, omegam, sigma8)
    jacsamples._weightsChanged()
    g.plot_2d([samples, s8samples, jacsamples], param_pair=pair, filled=True )
    g.export()


if True:
    g.plot_3d(samples, params=pair+['H0'], lims=ranges)

    #g.plot_3d('kids450fiducial', params=pair+['H0'], lims=ranges)

    g.add_2d_contours('base_' + s.defdata_all_lensing, param_pair=pair, filled=False) #, color=g.settings.solid_colors[1]

    g.add_2d_contours(s8samples, param_pair=pair, filled=False, color='red') #, color=g.settings.solid_colors[1]
#    g.add_2d_contours(s82, param_pair=pair, filled=False, color='green')  # , color=g.settings.solid_colors[1]

    #    g.add_2d_contours('fsigma-vel-theta', param_pair=pair, filled=False, color='yellow') #, color=g.settings.solid_colors[1]

    g.add_2d_contours('lensing_BAO', param_pair=pair, filled=False, color='orange', ls='--') #, color=g.settings.solid_colors[1]


    legends = [s.lensingall, r'DR12+ {$\theta$,$\Omega_b h^2$,$n_s$} prior', r'\textit{Planck} lensing+BAO+priors']
    g.add_legend(legends, legend_loc='upper right', colored_text=True, align_right=True)
    #gca().set_yticks(np.arange(0.75, 1, 0.05))
g.export()
