import numpy as np
from scipy.signal import fftconvolve
# convolve2D is just an alias to be consistent with convolve1D and allow playing with alternatives if desired
from scipy.signal import fftconvolve as convolve2D
from scipy import fftpack
from scipy.optimize import fsolve
import math

# numbers of the form 2^n3^m
fastFFT = np.array([ 2, 4, 8, 16, 32, 64, 96, 128, 144, 192, 256, 288, 384, 512, 576, 768, 1024, 1152, 1536, 2048, 2304, 3072,
4096 , 4608, 6144, 8192, 9216, 12288, 16384, 18432, 24576, 32768, 36864, 49152, 65536, 73728,
98304 , 131072, 147456, 196608, 262144, 294912, 393216, 524288, 589824, 786432, 1048576,
1179648 , 1572864, 2097152, 2359296, 3145728, 4194304, 4718592, 6291456, 8388608,
9437184 , 12582912, 16777216, 18874368, 25165824, 33554432, 37748736, 50331648, 67108864,
75497472 , 100663296, 134217728, 150994944, 201326592, 268435456, 301989888, 402653184, 452984832,
536870912 , 603979776, 805306368, 905969664], dtype=np.int)

rootpi = math.sqrt(np.pi)
root2pi = math.sqrt(2 * np.pi)
pisquared = np.pi ** 2

_kde_lmax_bandwidth = 7
_kde_consts = np.array([(1 + 0.5 ** (s + 0.5)) / 3 * np.prod(np.arange(1, 2 * s, 2)) / root2pi for s in  range(_kde_lmax_bandwidth - 1, 1, -1)])

def _bandwidth_fixed_point(t, N, I, logI, a2):
    # this implements the function t-zeta*gamma^[l](t)
    f = 2 * np.pi ** (2 * _kde_lmax_bandwidth) * np.dot(a2, np.exp(_kde_lmax_bandwidth * logI - I * (pisquared * t)))
    for s, const in zip(range(_kde_lmax_bandwidth - 1, 1, -1), _kde_consts):
        time = (2 * const / N / f) ** (2 / (3. + 2 * s));
        f = 2 * np.pi ** (2 * s) * np.dot(a2, np.exp(s * logI - I * (pisquared * time)));
    return t - (2 * N * rootpi * f) ** (-2. / 5);

def gaussian_kde_bandwidth(data, Neff, a=None):
    """
     Return optimal standard kernel bandwidth assuming data is binned from Neff independent samples
     Return value is the bandwidth in units of the range of data (i.e. multiply by max(data)-min(data))

     Uses Improved Sheather-Jones algorithm from
     Kernel density estimation via diffusion: Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
     Annals of Statistics, Volume 38, Number 5, pages 2916-2957. 
     http://arxiv.org/abs/1011.2602
    """
    I = np.arange(1, data.size) ** 2
    logI = np.log(I)
    if a is None: a = fftpack.dct(data / np.sum(data))
    a2 = (a[1:] / 2) ** 2
    t = 0.28 * Neff ** (-2. / 5)  # default value in case of failure
    try:
        t = fsolve(_bandwidth_fixed_point, t, (Neff, I, logI, a2), xtol=t / 20)
    except:
        print 'kde_bandwidth failed'
    return math.sqrt(t)


def convolve1D(x, y, mode):
    if min(x.shape[0], y.shape[0]) > 1000:
        return fftconvolve(x, y, mode)
    else:
        return np.convolve(x, y, mode)

def convolveGaussianDCT(x, sigma, pad_sigma=4, mode='same', cache={}):
    """
    1D convolution of x with Gaussian of width sigma pixels
    If pad_sigma>0, pads ends with zeros by int(pad_sigma*sigma) pixels
    Otherwise does unpadded fast cosine transform, hence reflection from the ends
    """
    global gauss, cache_s, cache_sigma

    fill = int(pad_sigma * sigma)
    actual_size = x.size + fill * 2
    if fill > 0:
        s = max(actual_size, fastFFT[np.searchsorted(fastFFT, actual_size)])
        fill2 = s - x.size - fill
        padded_x = np.pad(x, (fill, fill2), mode='constant')
    else:
        padded_x = x

    s = padded_x.size
    hnorm = sigma / float(s)
    gauss = cache.get((s, hnorm))
    if gauss is None:
        gauss = np.exp(-(np.arange(0, s) * (np.pi * hnorm)) ** 2 / 2.)
        cache[(s, hnorm)] = gauss
    res = fftpack.idct(fftpack.dct(padded_x, overwrite_x=fill > 0) * gauss, overwrite_x=fill > 0) / (2 * s)
    if fill == 0: return res
    if mode == 'same':
        return res[fill:-fill2]
    elif mode == 'valid':
        return res[fill * 2:-fill2 - fill]
    else: raise ValueError('mode not supported for convolveGaussianDCT')

def convolveGaussian(x, sigma, sigma_range=4, cache={}):
    """
    1D convolution of x with Gaussian of width sigma pixels
    x_max = int(sigma_range*sigma) the zero padding range at ends
    This uses periodic boundary conditions, and mode = 'same'
    """
    fill = int(sigma_range * sigma)
    actual_size = x.size + 2 * fill
    if fill > 0:
        s = max(actual_size, fastFFT[np.searchsorted(fastFFT, actual_size)])
    else:
        s = actual_size
    gauss = cache.get((fill, actual_size, sigma))
    if gauss is None:
        hnorm = sigma / float(s)
        ps = np.arange(1, s + 1) // 2
        gauss = np.exp(-(ps * (np.pi * hnorm)) ** 2 * 2)
        cache[(fill, actual_size, sigma)] = gauss
    res = fftpack.irfft(fftpack.rfft(x, s) * gauss, s)
    return res[:x.size]


def convolveGaussianTrunc(x, sigma, sigma_range=4, mode='same', cache={}):
    """
    1D convolution of x with Gaussian of width sigma pixels
    x_max = int(sigma_range*sigma) determines the finite support (in pixels) of the truncated gaussian
    This uses normalized finite range approximation to Gaussian
    """
    fill = int(sigma_range * sigma)
    actual_size = x.size + 2 * fill
    if fill > 0:
        s = max(actual_size, fastFFT[np.searchsorted(fastFFT, actual_size)])
    else:
        s = actual_size
    gauss = cache.get((fill, actual_size, sigma))
    if gauss is None:
        points = np.arange(-fill, fill + 1)
        Win = np.exp(-(points / sigma) ** 2 / 2.)
        Win /= np.sum(Win)
        gauss = np.fft.rfft(Win, s)
        cache[(fill, actual_size, sigma)] = gauss
    res = np.fft.irfft(np.fft.rfft(x, s) * gauss, s)[:actual_size]
    if mode == 'same':
        return res[fill:-fill]
    elif mode == 'full':
        return res
    elif mode == 'valid':
        return res[2 * fill:-2 * fill]
