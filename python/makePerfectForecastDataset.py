# take CAMB file (e.g. test_lensedCls.dat') and produce dataset with given noise for testing
from numpy import *
import iniFile
import shutil
import sys, os



# # Edit parameters you want to change here

lensedTotClFileRoot = 'flatLCDM_lensedCls'
reconNoise = 'reconNoise.dat'
lensPotential = 'flatLCDM_lenspotentialCls.dat'
fwhm_arcmin = 7.
# Noise var is N_l in muK^2 for white noise
NoiseVar = 2e-4
ENoiseFac = 2
lmin = 2
lmax = 2000
fsky = 1

DoLensing = True


outDir = 'data/'
# os.path.dirname(sys.path[0])+'/data/'
ini = iniFile.iniFile()
dataset = ini.params

# change this if you don't want to use all pol
dataset['fields_use'] = 'T E B'


# #Now produce the Planck_like files
fwhm = fwhm_arcmin / 60
xlc = 180 * sqrt(8.*log(2.)) / pi
sigma2 = (fwhm / xlc) ** 2

outPath = '%DATASETDIR%'
outRoot = lensedTotClFileRoot + '_exactsim'
NoiseOut = []
for l in range(lmax + 1):
    NoiseCl = l * (l + 1) / 2 / pi * NoiseVar * exp(l * (l + 1) * sigma2)
    NoiseOut.append([NoiseCl, ENoiseFac * NoiseCl, ENoiseFac * NoiseCl])


outfile = open(outDir + outRoot + '_Noise.dat', 'w')
for l in range(2, lmax + 1):
    outfile.write("%d " % l + " ".join("%E" % elem for elem in NoiseOut[l]) + "\n")
outfile.close()

dataset['fullsky_exact_fksy'] = fsky
dataset['name'] = outRoot
dataset['dataset_format'] = 'CMBLike'
dataset['lowl_exact'] = False
dataset['highl_cl'] = True
dataset['like_approx'] = 3


dataset['bin_width'] = 1
dataset['cl_lmin'] = lmin
dataset['cl_lmax'] = lmax

dataset['cl_hat_includes_noise'] = False
dataset['cl_hat_includes_pointsources'] = False
dataset['cl_hat_file'] = outPath + outRoot + '.dat'
dataset['cl_hat_order'] = 'TT EE BB TE'
dataset['cl_noise_file '] = outPath + outRoot + '_Noise.dat'
dataset['cl_noise_order'] = 'TT EE BB'
dataset['cl_noise_includes_pointsources'] = False
dataset['point_source_cl'] = ''
dataset['beam_modes_number'] = 0

if DoLensing:
    dataset['lensing_recon_ncl'] = 1
    dataset['cl_hat_phi_file'] = outPath + outRoot + '_PhiEst.dat'
    dataset['cl_noise_phi_file'] = outPath + outRoot + '_PhiNoise.dat'
    dataset['cl_hat_includes_noise'] = False
    shutil.copy(reconNoise, outDir + outRoot + '_PhiNoise.dat')
    fout = open(outDir + outRoot + '_PhiEst.dat', 'w')
    f = open(lensPotential)
    for line in f:
            row = line.split()
            fout.write(str(int(row[0])) + ' ' + str(row[5]) + '\n')
    f.close()
    fout.close()


ini.saveFile(outDir + outRoot + '.dataset')

shutil.copy(lensedTotClFileRoot + '.dat', outDir + outRoot + '.dat')
