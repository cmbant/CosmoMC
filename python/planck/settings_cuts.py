from __future__ import absolute_import
import re
import copy

from paramgrid import batchjob
from six.moves import range


# Look at Alens=1 for fixed lensing template amplitude

ini_dir = 'batch2/'

defaults = ['common.ini']

importanceDefaults = ['importance_sampling.ini']

override_defaults = ['pico.ini']

Camspec = 'CAMspec_defaults.ini'
highL = 'highL'
lowl = 'lowl'
# dataset names
tauprior = {'prior[tau]': '0.07 0.015'}
tauname = 'tau07'
WMAPtau = {'prior[tau]': '0.09 0.013'}

planck_detsets = ['nonclik_detsets.ini']
planck_CS = ['nonclik.ini']

CamSpec = [
    batchjob.dataSet(['CamSpecHM', 'TT', 'lowl', tauname], ['nonclik.ini', 'CAMspec_TT.ini', 'lowl.ini', tauprior])]

start_at_bestfit = False
newCovmats = True

lmaxs = list(range(700, 2501, 50))

groups = []

g = batchjob.jobGroup('main')
# Main group with just tau prior
g.datasets = []
for lmax in lmaxs:
    for d in copy.deepcopy(CamSpec):
        d.add('lmax' + str(lmax), {'camspec_lmax': (str(lmax) + ' ') * 6})
        g.datasets.append(d)
g.params = [[], ['Alens']]
groups.append(g)

covrenames = []
covrenames.append(['_tau07_lowl', '_lowTEB'])
covrenames.append(['_tau07', '_lowTEB'])

covNameMappings = {'v910CMH': 'CamSpec', 'CamSpecHM': 'CamSpec', 'plikHMv17': 'plik', 'plikDSv16sz': 'plik',
                   'plikHMv16sz': 'plik', 'plikHMv16bin1sz': 'plik', 'bin1l80sz': 'plik', 'tau07': 'lowTEB',
                   'no217auto': '', 'no217': '', 'no143': ''}

for mx in lmaxs:
    covNameMappings['lmax' + str(mx)] = ''


def covRenamer(name):
    # renamed = re.sub(r'_lmax.*', '', name, re.I)
    renamed = re.sub(r'_v.*', '_CamSpecHM_TT_lowTEB', name, re.I)

    if renamed == name:
        return []
    else:
        return [renamed]
