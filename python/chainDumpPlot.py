from __future__ import absolute_import
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import numpy as np
import os
import fnmatch
from getdist import types

"""
    load and plot chain cl data, as export by cosmomc using action =1 and redo_output_txt_theory=T, redo_output_txt_dir=xxx
"""

colorbar_label_rotation_angle = -90
colorbar_tick_label_vertical = False
label_dict = {}
color_ticks = None
color_ticklabels = None


def loadFiles(pat, max_files=400, cl_ix=1, calParam=None, rescale=1):
    d, name = os.path.split(pat + '_?_*')
    pars = []
    cls = []
    ls = None
    for f in os.listdir(d):
        if fnmatch.fnmatch(f, name + '.theory_cl'):
            fname = f.replace('.theory_cl', '')
            params = types.BestFit(os.path.join(d, fname + '.pars'))
            pars.append(params)
            cl = np.loadtxt(os.path.join(d, fname + '.theory_cl'))
            if calParam:
                cal = params.parWithName(calParam).best_fit
            else:
                cal = 1
            cls.append(cl[:, cl_ix] / cal ** 2 * rescale)
            if ls is None: ls = cl[:, 0]
            if len(pars) == max_files: break
    if not len(pars): raise Exception('No files found: ' + pat)
    return ls, pars, cls


def CLdensity(samples, param, scalarMap, nbins=300, deltaL=1, Lmax=200, color_pow=0.25):
#    H, xedges, yedges = np.histogram2d(x, y, bins=(100, 20), weights=param)
    minCl = np.min(samples) / 1.04
    if minCl > 0: minCl = 0
    maxCl = np.max(samples) * 1.04
    delta = (maxCl - minCl) / (nbins - 1)
    nL = (Lmax - 1) / deltaL
    sums = np.zeros((nL, nbins))
    psums = np.zeros((nL, nbins))

    for i in  range(param.shape[0]):
        indices = np.floor((samples[i, ::deltaL] - minCl) / delta).astype(int)
        for j in range(nL):
            psums[j, indices[j]] += param[i]
            sums[j, indices[j]] += 1

    window = np.arange(-2.5, 2.5 + 0.1, 0.5)
    window = np.exp(-abs(window ** 3) / 2)
    window /= np.sum(window)
    for i in range(nL):
        counts = sums[i, :]
        psums[i, :] = np.convolve(psums[i, :], window, 'same')
        sums[i, :] = np.convolve(counts, window, 'same')
        positive = np.where(counts > 0)
        psums[i, positive] /= counts[positive]
        sums[i, :] /= np.max(sums[i, :])

        maxix = np.max(positive)
        minix = np.min(positive)
        psums[i, maxix:] = psums[i, maxix]
        psums[i, :minix] = psums[i, minix]

    if deltaL == 1:
        for i in range(nbins):
            sums[:, i] = np.convolve(sums[:, i] , window, 'same')

    colorVal = np.zeros((nbins, nL, 4))
    for i in range(nL):
            for j in range(nbins):
                colorVal[j, i, :] = scalarMap.to_rgba(psums[i, j], sums[i, j] ** color_pow)
    im = plt.imshow(colorVal, origin='lower', aspect='auto', extent=[2, Lmax, 0, maxCl])
    im.set_interpolation('bicubic')


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=200):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def getColorMap(parvals):
    jet = plt.get_cmap('jet')
    cmap = truncate_colormap(jet, 0.05, 0.9)
    # delta = np.max(parvals) - np.min(parvals)
    delta = 0
    cNorm = colors.Normalize(vmin=np.min(parvals) + delta / 10, vmax=np.max(parvals) - delta / 10)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    scalarMap.set_array(parvals)
    return scalarMap, cNorm

def setLabels(scalarMap, cNorm, par):
    plt.xlabel('$L$')
    cb = plt.colorbar(scalarMap, norm=cNorm)
    lab = label_dict.get(par.name, par.label)

    cb.set_label('$' + lab + '$', rotation=colorbar_label_rotation_angle, labelpad=10)
    if colorbar_tick_label_vertical:
        if color_ticks is not None:
            cb.set_ticks(color_ticks)
        if color_ticklabels is not None:
            cb.set_ticklabels(color_ticklabels)
        for ticklabel in cb.ax.get_yticklabels():
            ticklabel.set_rotation(-90)
        if color_ticks is None:
            labels = [label.get_text() for label in cb.ax.yaxis.get_ticklabels()[::2]]
            cb.ax.yaxis.set_ticks(cb.ax.yaxis.get_ticklocs()[::2])
            cb.ax.yaxis.set_ticklabels(labels)


def rainbowLinePlot(ls, pars, cls, parname, lpow=2, delta_to=None, last_colors=None, alpha=1):

    parvals = [par.parWithName(parname).best_fit for par in pars]
    scalarMap, cNorm = last_colors or getColorMap(parvals)
    sc = (ls * (ls + 1.)) ** (lpow / 2. - 1)
    for parval, cl in zip(parvals, cls):
        colorVal = scalarMap.to_rgba(parval)
        if delta_to is not None:
            plt.plot(ls, (cl - delta_to) * sc, color=colorVal, alpha=alpha)
        else:
            plt.plot(ls, cl * sc, color=colorVal, alpha=alpha)
    setLabels(scalarMap, cNorm, pars[0].parWithName(parname))
    return scalarMap, cNorm


def rainbowDensityPlot(ls, pars, cls, parname, lpow=2, delta_to=None, Lmax=200, color_pow=0.25, last_colors=None):

    parvals = np.array([par.parWithName(parname).best_fit for par in pars])
    scalarMap, cNorm = last_colors or getColorMap(parvals)
    samples = np.zeros((len(cls), Lmax - 1))
    for i, cl in enumerate(cls):
        samples[i, :] = cl[:Lmax - 1]
    CLdensity(samples, parvals, scalarMap, Lmax=Lmax, color_pow=color_pow)
    setLabels(scalarMap, cNorm, pars[0].parWithName(parname))
    return scalarMap, cNorm
