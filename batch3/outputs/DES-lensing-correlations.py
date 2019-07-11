# Plot DES lensing bandpowers. Need planck bestfit minimum run separately

import planckStyle as s

import numpy as np
import pylab as plt
import sys
import tempfile
import hashlib

sys.path.insert(0, r'C:\Work\Dist\git\camb')
from cosmomc_to_camb import get_camb_params
import camb
from planck.DES import DES_like
from getdist.types import BestFit
import os
import pickle
from copy import deepcopy

g = s.getSubplotPlotter()
fast_test = False
# camb.set_halofit_version('mead')

DESdataset = os.path.join(os.path.dirname(__file__), '../../data/DES/DES_1YR_final.dataset')

figsize = 6
plots = ['lens']


# plots=['auto','gammt','lens']

def getdic(paramdic):
    paramdic['tau'] = 0.055
    paramdic['hslofit_version'] = 'mead'
    return paramdic


if 'lens' in plots:
    DESlens = DES_like(DESdataset, dataset_params={'used_data_types': 'xip xim'})

    samples = g.sampleAnalyser.samplesForRoot('base_DESlens_DESpriors')
    planck_fit_lens = getdic(BestFit(r'C:\Tmp\Planck\2017\DES-Planck_bf\base_DESlens_DESpriors_planckbf.minimum',
                                     want_fixed=True).getParamDict())
    lens_theory = DESlens.get_theory_for_params(planck_fit_lens)

    planck_fit_lens['DES_AIA'] = 0
    lens_theory_noIA = DESlens.get_theory_for_params(planck_fit_lens)
    hash = hashlib.md5(('%s%s%s'%(fast_test, samples.jobItem.chainRoot,
                        os.path.getmtime(samples.jobItem.chainRoot + '_1.txt'))).encode('ascii')).hexdigest()
    cachename = os.path.join(tempfile.gettempdir(), 'planck2018' + hash + '.DESlens_samps')
    if os.path.isfile(cachename):
        print('reading cache %s' % cachename)
        with open(cachename, 'rb') as inp:
            theory_samps = pickle.load(inp)
    else:
        print('Calculating DES spectrum samples')
        theory_samps = [np.zeros((DESlens.nzbins, DESlens.nzbins), dtype=object),
                        np.zeros((DESlens.nzbins, DESlens.nzbins), dtype=object)]

        for f1, f2 in DESlens.bin_pairs[0]:
            theory_samps[0][f1, f2] = []
            theory_samps[1][f1, f2] = []

        if fast_test:
            paramdic = getdic(samples.getParamSampleDict(0))
            camb_pars = get_camb_params(paramdic)
            results, PKdelta, PK = DESlens.get_camb_theory(camb_pars)
        ixs = samples.randomSingleSamples_indices()[::4]
        for ix, i in enumerate(ixs):
            print("%s of %s" % (ix, len(ixs)))
            paramdic = getdic(samples.getParamSampleDict(i))
            if fast_test:
                th = np.array(
                    DESlens.get_theory_for_params(paramdic, camb_pars=camb_pars, camb_results=[results, PKdelta, PK]))
            else:
                th = np.array(DESlens.get_theory_for_params(paramdic))
            for f1, f2 in DESlens.bin_pairs[0]:
                theory_samps[0][f1, f2].append(th[0][f1, f2])
                theory_samps[1][f1, f2].append(th[1][f1, f2])
        with open(cachename, 'wb') as output:
            print('cacheing %s' % cachename)
            pickle.dump(theory_samps, output, pickle.HIGHEST_PROTOCOL)

    mean_th = deepcopy(lens_theory)
    lims = deepcopy(lens_theory)

    for f1, f2 in DESlens.bin_pairs[0]:
        arr = np.array(theory_samps[0][f1, f2])
        mean_th[0][f1, f2] = np.mean(arr, axis=0)
        lims[0][f1, f2] = np.percentile(arr, [2.5, 16, 84, 97.5], axis=0)
        arr = np.array(theory_samps[1][f1, f2])
        mean_th[1][f1, f2] = np.mean(arr, axis=0)
        lims[1][f1, f2] = np.percentile(arr, [2.5, 16, 84, 97.5], axis=0)

    g.fig, axs = plt.subplots(DESlens.nzbins, DESlens.nzbins, figsize=(figsize, figsize), sharex='col', sharey='row')
    for ax in axs.reshape(-1): ax.axis('off')
    for f1, f2 in DESlens.bin_pairs[0]:
        ax = axs[f2, f1]
        ax.axis('on')
        xip = DESlens.data_arrays[0][f1, f2]
        fac = 1e4 * DESlens.theta_bins
        xip = xip * fac
        ax.errorbar(DESlens.theta_bins, xip, fac * DESlens.errors[0][f1, f2], color='C0', fmt='.')
        ax.semilogx(DESlens.theta_bins, fac * lens_theory[0][f1, f2], color='k', ls='--')
        ax.semilogx(DESlens.theta_bins, fac * (lens_theory[0][f1, f2] - lens_theory_noIA[0][f1, f2]), color='k', ls=':')

        ax.fill_between(DESlens.theta_bins, fac * lims[0][f1, f2][0, :], fac * lims[0][f1, f2][3, :], color='C2',
                        alpha=0.2)
        ax.fill_between(DESlens.theta_bins, fac * lims[0][f1, f2][1, :], fac * lims[0][f1, f2][2, :], color='C2',
                        alpha=0.5)

        ax.axhline(0, color='k', lw=0.5)
        ax.text(0.1, 0.85, '%s-%s' % (f1 + 1, f2 + 1), transform=ax.transAxes,
                bbox=dict(facecolor='ivory', alpha=0.5, boxstyle="square,pad=0.25"), fontsize=7)
        ax.axvspan(DESlens.theta_bins[0] - 0.5, DESlens.ranges['xip'][f1][f2][0], color='gray', alpha=0.1, zorder=-100)
        ax.set_xlim(DESlens.theta_bins[0] - 0.5, None)
        ax.get_xaxis().set_tick_params(which='both', direction='in', right=False, top=False)
        ax.get_yaxis().set_tick_params(which='both', direction='in', right=False, top=False)

    plt.subplots_adjust(wspace=0, hspace=0)
    g.fig.text(0.06, 0.5, r'$10^4 \theta\xi_+(\theta)\, [{\rm arcmin}]$', va='center', rotation='vertical', fontsize=12)
    g.fig.text(0.5, 0.06, r'$\theta\, [{\rm arcmin}]$', ha='center', fontsize=12)
    g.fig.savefig('../../outputs/DES-shear-spec-plus.pdf', bbox_inches='tight')

    g.fig, axs = plt.subplots(DESlens.nzbins, DESlens.nzbins, figsize=(figsize, figsize), sharex='col', sharey='row')
    for ax in axs.reshape(-1): ax.axis('off')
    for f1, f2 in DESlens.bin_pairs[0]:
        ax = axs[f2, f1]
        ax.axis('on')
        xim = DESlens.data_arrays[1][f1, f2]
        fac = 1e4 * DESlens.theta_bins
        facm = fac
        xim = xim * fac
        ax.errorbar(DESlens.theta_bins, xim, facm * DESlens.errors[1][f1, f2], color='C0', fmt='.')
        ax.semilogx(DESlens.theta_bins, facm * lens_theory[1][f1, f2], color='k', ls='--')
        ax.semilogx(DESlens.theta_bins, fac * (lens_theory[1][f1, f2] - lens_theory_noIA[1][f1, f2]), color='k', ls=':')

        ax.fill_between(DESlens.theta_bins, fac * lims[1][f1, f2][0, :], fac * lims[1][f1, f2][3, :], color='C2',
                        alpha=0.2)
        ax.fill_between(DESlens.theta_bins, fac * lims[1][f1, f2][1, :], fac * lims[1][f1, f2][2, :], color='C2',
                        alpha=0.5)
        ax.axhline(0, color='k', lw=0.5)
        ax.text(0.1, 0.85, '%s-%s' % (f1 + 1, f2 + 1), transform=ax.transAxes,
                bbox=dict(facecolor='ivory', alpha=0.5, boxstyle="square,pad=0.25"), fontsize=7)
        ax.axvspan(DESlens.theta_bins[0] - 0.5, DESlens.ranges['xim'][f1][f2][0], color='gray', alpha=0.1, zorder=-100)
        ax.set_xlim(DESlens.theta_bins[0] - 0.5, None)
        ax.get_xaxis().set_tick_params(which='both', direction='in', right=False, top=False)
        ax.get_yaxis().set_tick_params(which='both', direction='in', right=False, top=False)

    plt.subplots_adjust(wspace=0, hspace=0)
    g.fig.text(0.06, 0.5, r'$10^4 \theta\xi_-(\theta)\, [{\rm arcmin}]$', va='center', rotation='vertical',
               fontsize=12)
    g.fig.text(0.5, 0.06, r'$\theta\, [{\rm arcmin}]$', ha='center', fontsize=12)
    g.fig.savefig('../../outputs/DES-shear-spec-minus.pdf', bbox_inches='tight')

if 'auto' in plots or 'gammat' in plots:

    DES = DES_like(DESdataset)

    samples = g.sampleAnalyser.samplesForRoot('base_DES_DESpriors')
    planck_fit = getdic(BestFit(r'C:\Tmp\Planck\2017\DES-Planck_bf\base_DESwt_planckbf.minimum',
                                want_fixed=True).getParamDict())
    theory = DES.get_theory_for_params(planck_fit)

    planck_fit['DES_AIA'] = 0
    theory_noIA = DES.get_theory_for_params(planck_fit)

    cachename = os.path.join(tempfile.gettempdir(), 'planck2018' +
                             str(hash((fast_test, samples.jobItem.chainRoot,
                                       os.path.getmtime(samples.jobItem.chainRoot + '_1.txt'))))
                             + '.DES_samps')
    if os.path.isfile(cachename):
        print('reading cache %s' % cachename)
        with open(cachename, 'rb') as inp:
            auto_theory_samps, galgal_theory_samps = pickle.load(inp)
    else:
        print('Calculating DES spectrum samples')
        auto_theory_samps = np.zeros((DES.nwbins, DES.nwbins), dtype=object)
        for f1, f2 in DES.bin_pairs[3]:
            auto_theory_samps[f1, f2] = []
        galgal_theory_samps = np.zeros((DES.nwbins, DES.nzbins), dtype=object)
        for f1, f2 in DES.bin_pairs[2]:
            galgal_theory_samps[f1, f2] = []

        if fast_test:
            paramdic = getdic(samples.getParamSampleDict(0))
            camb_pars = get_camb_params(paramdic)
            results, PKdelta, PK = DES.get_camb_theory(camb_pars)
        ixs = samples.randomSingleSamples_indices()[::4]
        for ix, i in enumerate(ixs):
            print("%s of %s" % (ix, len(ixs)))
            paramdic = getdic(samples.getParamSampleDict(i))
            if fast_test:
                th = np.array(
                    DES.get_theory_for_params(paramdic, camb_pars=camb_pars, camb_results=[results, PKdelta, PK]))
            else:
                th = np.array(DES.get_theory_for_params(paramdic))
            for f1, f2 in DES.bin_pairs[2]:
                galgal_theory_samps[f1, f2].append(th[2][f1, f2])
            for f1, f2 in DES.bin_pairs[3]:
                auto_theory_samps[f1, f2].append(th[3][f1, f2])

        with open(cachename, 'wb') as output:
            print('cacheing %s' % cachename)
            pickle.dump([auto_theory_samps, galgal_theory_samps], output, pickle.HIGHEST_PROTOCOL)

    mean_th = deepcopy(theory)
    lims = deepcopy(theory)

    if 'gammat' in plots:
        # galaxy galaxy cross
        for f1, f2 in DES.bin_pairs[2]:
            arr = np.array(galgal_theory_samps[f1, f2])
            mean_th[2][f1, f2] = np.mean(arr, axis=0)
            lims[2][f1, f2] = np.percentile(arr, [2.5, 16, 84, 97.5], axis=0)

        g.fig, axs = plt.subplots(DES.nzbins, DES.nwbins, figsize=(figsize, figsize * 4. / 5), sharex='col',
                                  sharey='row')
        for f1, f2 in DES.bin_pairs[2]:
            ax = axs[f2, f1]
            corr = DES.data_arrays[2][f1, f2]
            fac = 1e2 * DES.theta_bins
            corr *= fac
            ax.errorbar(DES.theta_bins, corr, fac * DES.errors[2][f1, f2], color='C0', fmt='.')
            ax.semilogx(DES.theta_bins, fac * theory[2][f1, f2], color='k', ls='--')
            ax.semilogx(DES.theta_bins, fac * (theory[2][f1, f2] - theory_noIA[2][f1, f2]), color='k', ls=':')

            ax.fill_between(DES.theta_bins, fac * lims[2][f1, f2][0, :], fac * lims[2][f1, f2][3, :], color='C2',
                            alpha=0.2)
            ax.fill_between(DES.theta_bins, fac * lims[2][f1, f2][1, :], fac * lims[2][f1, f2][2, :], color='C2',
                            alpha=0.5)

            ax.axhline(0, color='k', lw=0.5)
            ax.text(0.1, 0.85, '%s-%s' % (f1 + 1, f2 + 1), transform=ax.transAxes,
                    bbox=dict(facecolor='ivory', alpha=0.5, boxstyle="square,pad=0.25"), fontsize=7)
            ax.axvspan(DES.theta_bins[2] - 0.5, DES.ranges['gammat'][f1][f2][0], color='gray', alpha=0.1, zorder=-100)
            ax.set_xlim(DES.theta_bins[2] - 0.5, None)
            ax.get_xaxis().set_tick_params(which='both', direction='in', right=False, top=False)
            ax.get_yaxis().set_tick_params(which='both', direction='in', right=False, top=False)

        plt.subplots_adjust(wspace=0, hspace=0)
        g.fig.text(0.06, 0.5, r'$10^2 \gamma_t(\theta)\, [{\rm arcmin}]$', va='center', rotation='vertical',
                   fontsize=12)
        g.fig.text(0.5, 0.04, r'$\theta\, [{\rm arcmin}]$', ha='center', fontsize=12)
        g.fig.savefig('../../outputs/DES-galgal-spec.pdf', bbox_inches='tight')

    if 'auto' in plots:
        # galaxy autocorrelation
        DESauto = DES_like(DESdataset, dataset_params={'used_data_types': 'wtheta'})

        planck_fit = getdic(BestFit(r'C:\Tmp\Planck\2017\DES-Planck_bf\base_DESauto_planckbf.minimum',
                                    want_fixed=True).getParamDict())

        theory = DES.get_theory_for_params(planck_fit)

        for f1 in range(DES.nwbins):
            arr = np.array(auto_theory_samps[f1, f1])
            mean_th[3][f1, f1] = np.mean(arr, axis=0)
            lims[3][f1, f1] = np.percentile(arr, [2.5, 16, 84, 97.5], axis=0)
        g.fig, axs = plt.subplots(1, DES.nwbins, figsize=(figsize, figsize / 4.5), sharex='col', sharey='row')
        for f1, ax in zip(range(DES.nwbins), axs.reshape(-1)):
            corr = DES.data_arrays[3][f1, f1]
            fac = DES.theta_bins
            corr *= fac
            ax.errorbar(DES.theta_bins, corr, fac * DES.errors[3][f1, f1], color='C0', fmt='.')
            ax.semilogx(DES.theta_bins, fac * theory[3][f1, f1], color='k', ls='--')

            ax.fill_between(DES.theta_bins, fac * lims[3][f1, f1][0, :], fac * lims[3][f1, f1][3, :], color='C2',
                            alpha=0.2)
            ax.fill_between(DES.theta_bins, fac * lims[3][f1, f1][1, :], fac * lims[3][f1, f1][2, :], color='C2',
                            alpha=0.5)

            ax.axhline(0, color='k', lw=0.5)
            ax.text(0.1, 0.85, '%s-%s' % (f1 + 1, f1 + 1), transform=ax.transAxes,
                    bbox=dict(facecolor='ivory', alpha=0.5, boxstyle="square,pad=0.25"), fontsize=7)
            ax.axvspan(DES.theta_bins[3] - 0.5, DES.ranges['wtheta'][f1][f1][0], color='gray', alpha=0.1, zorder=-100)
            ax.set_xlim(DES.theta_bins[3] - 0.5, None)
            ax.get_xaxis().set_tick_params(which='both', direction='in', right=False, top=False)
            ax.get_yaxis().set_tick_params(which='both', direction='in', right=False, top=False)
            if f1 == 0:
                ax.set_ylabel(r'$\theta w(\theta)\, [{\rm arcmin}]$', fontsize=12)
            if f1 == 2:
                ax.set_xlabel(r'$\theta\, [{\rm arcmin}]$', fontsize=12)

        plt.subplots_adjust(wspace=0, hspace=0)
        g.fig.savefig('../../outputs/DES-auto-spec.pdf', bbox_inches='tight')
