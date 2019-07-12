# Get python camb parameters from dictionary of cosmomc parameter values
from __future__ import absolute_import
from __future__ import print_function
import os
import sys

try:
    import camb
    from camb import model, initialpower

    print('Using CAMB at: %s' % os.path.dirname(camb.__file__))
except ImportError as e:
    print('error importing CAMB, run setup.py install first', e)


def get_camb_params(p, num_massive_neutrinos=1, neutrino_hierarchy='degenerate', inpars=None, halofit_version = 'mead'):
    pars = inpars or camb.CAMBparams()
    if p.get('wa', 0) or p.get('alpha1', 0) or p.get('Alens', 1) != 1 or p.get('Aphiphi', 1) != 1:
        raise Exception('Parameter not currrently supported by Cosmomc-> camb converter')
    pars.set_cosmology(H0=p['H0'], ombh2=p['omegabh2'], omch2=p['omegach2'], mnu=p.get('mnu', 0.06),
                       omk=p.get('omegak', 0),
                       tau=p['tau'], deltazrei=p.get('deltazrei', None), nnu=p.get('nnu', 3.046),
                       YHe=p.get('yheused', None),
                       meffsterile=p.get('meffsterile', 0),
                       num_massive_neutrinos=num_massive_neutrinos, neutrino_hierarchy=neutrino_hierarchy)
    pars.InitPower.set_params(ns=p['ns'], r=p.get('r', 0), As=p['A'] * 1e-9, nrun=p.get('nrun', 0),
                              nrunrun=p.get('nrunrun', 0))
    pars.set_dark_energy(w=p.get('w', -1))
    pars.NonLinearModel.set_params(halofit_version = halofit_version)
    pars.set_for_lmax(2500, lens_potential_accuracy=1)
    return pars
