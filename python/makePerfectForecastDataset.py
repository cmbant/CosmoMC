# take CAMB file (e.g. test_lensedCls.dat) and produce dataset with given noise for testing
# Use in cosmomc .ini file using e.g.
# cmb_dataset[MyForecast]=data/MyForecast/test_lensedCls_exactsim.dataset

from __future__ import absolute_import
import shutil
import os
import numpy as np
from getdist import IniFile

# Edit parameters you want to change here

root = r'test_lensedCls'

lensedTotClFileRoot = os.path.join(os.path.dirname(__file__), '..', 'data', 'base_plikHM_TT_lowTEB.minimum.theory_cl')
dataCols = 'TT TE EE BB PP'
outDir = 'z:\\'

# these numbers are Planck-like
# Noise var is N_l in muK^2 for white noise
# note  NoiseVar = (muKArcmin * np.pi / 180 / 60.) ** 2

fwhm_arcmin = 5.
# NoiseVar = 2e-4
NoiseVar = 4e-5
# Pol noise var = ENoiseFac * NoiseVar
# 2 normally, but for Planck only half detectors are polarized
ENoiseFac = 4
lmin = 2
lmax = 2500
fsky = 0.57


# os.path.dirname(sys.path[0])+'/data/'
ini = IniFile()
dataset = ini.params

# change this if you don't want to use all pol
dataset['fields_use'] = 'T E'


# #Now produce the Planck_like files
fwhm = fwhm_arcmin / 60
xlc = 180 * np.sqrt(8.*np.log(2.)) / np.pi
sigma2 = (fwhm / xlc) ** 2

outPath = ''
outRoot = root + '_exactsim'
NoiseOut = []

# d = loadtxt(lensedTotClFileRoot)
# SN = 0
for l in range(lmax + 1):
    NoiseCl = l * (l + 1.) / 2 / np.pi * NoiseVar * np.exp(l * (l + 1) * sigma2)
    NoiseOut.append([NoiseCl, ENoiseFac * NoiseCl, ENoiseFac * NoiseCl])
#    if (l >= 2): SN += (2 * l + 1) * fsky * (d[l - 2, 1] / (NoiseCl + d[l - 2, 1])) ** 2
# print 'Number of modes: ', SN

outfile = open(outDir + outRoot + '_Noise.dat', 'w')
for l in range(2, lmax + 1):
    outfile.write("%d " % l + " ".join("%E" % elem for elem in NoiseOut[l]) + "\n")
outfile.close()

dataset['fullsky_exact_fksy'] = fsky
dataset['name'] = outRoot
dataset['dataset_format'] = 'CMBLike2'
dataset['like_approx'] = 'exact'

dataset['cl_lmin'] = lmin
dataset['cl_lmax'] = lmax

dataset['binned'] = False


dataset['cl_hat_includes_noise'] = False

dataset['cl_hat_file'] = outPath + outRoot + '.dat'
dataset['cl_hat_order'] = dataCols
dataset['cl_noise_file '] = outPath + outRoot + '_Noise.dat'
dataset['cl_noise_order'] = 'TT EE BB'


ini.saveFile(outDir + outRoot + '.dataset')

shutil.copy(lensedTotClFileRoot, outDir + outRoot + '.dat')
