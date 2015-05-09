from __future__ import absolute_import
from __future__ import print_function
import os

from . import plik_postprocess
from getdist import inifile
from paramgrid import batchjob_args


Opts = batchjob_args.batchArgs('add plik params and bestfits', importance=True)

Opts.parser.add_argument('--finished', action='store_true', help='only run on completed chains')

(batch, args) = Opts.parseForBatch()


for jobItem in Opts.filteredBatchItems():
    if args.finished and not jobItem.chainFinished(): continue
    name = jobItem.chainRoot + '.paramnames'
    properties = jobItem.propertiesIni()
    if 'plik' in name and os.path.exists(name) and not properties.bool('plik_foregrounds', False):

        ini = inifile.IniFile(jobItem.chainRoot + '.inputparams')
        dat = ini.string('clik_data_plik', '')
        params = ini.string('clik_params_plik', '')
        hasderived = dat
        if not dat:
            dat = ini.string('clik_data_plikTE', '')
            params = ini.string('clik_params_plikTE', '')
        if not dat:
            dat = ini.string('clik_data_plikEE', '')
            params = ini.string('clik_params_plikEE', '')
        if not dat: raise Exception('no clik file found:' + jobItem.chainRoot)
        dat = dat.replace('%DATASETDIR%', './data/')
        params = params.replace('%DATASETDIR%', './data/')
        print(dat, params)

        minroot = jobItem.chainRoot + '.minimum'
        if os.path.exists(jobItem.chainRoot + '.minimum'):
            print(minroot)
            plik_postprocess.main_fgfile([ '', dat, params, minroot, minroot + '.plik_foregrounds'])
        else: minroot = jobItem.chainRoot + '.ranges'
        chains = jobItem.chainNames()
        if hasderived: plik_postprocess.main_fg2000([ '', dat, params, minroot, jobItem.chainRoot + '.paramnames'] + chains)
        properties.params['plik_foregrounds'] = True
        properties.saveFile()

print('Done. Re-run getdist to update getdist outputs.')
