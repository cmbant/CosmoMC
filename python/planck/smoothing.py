import numpy as np
import os
from scipy import signal
import numpy as np

def_smooth = 40


def gaussian_smooth_1d(data, smooth=def_smooth, nsigma=5):
    win = signal.gaussian(nsigma * smooth, std=smooth)
    res = np.convolve(data, win, mode='same')
    return res / np.convolve(np.ones(data.shape), win, mode='same')


def smoothed_cov(cov, smooth=def_smooth):
    scov = cov.copy()
    for i in range(cov.shape[0]):
        scov[i, :] = gaussian_smooth_1d(cov[i, :], smooth)
    scov2 = scov.copy()
    for i in range(cov.shape[0]):
        scov2[:, i] = gaussian_smooth_1d(scov[:, i], smooth)
    return scov2
