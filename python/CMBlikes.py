# Load CosmoMC format .dataset files with lensing likelihood data
# AL July 2014
# note this is not well tested with final published versions of likelihoods
# Does not handle calibration parameter
from __future__ import absolute_import
from __future__ import print_function
from matplotlib import pyplot as plt
import os
import numpy as np
import sys
from getdist import IniFile


def readTextCommentColumns(fname, cols):
        with open(fname) as f:
            x = f.readline().strip()
            if x[0] != '#': raise Exception('No Comment')
        incols = x[1:].split()
        colnums = [incols.index(col) for col in cols]
        return np.loadtxt(fname, usecols=colnums, unpack=True)

def readWithHeader(fname):
        with open(fname) as f:
            x = f.readline().strip()
            if x[0] != '#': raise Exception('No Comment')
            x = x[1:].split()
        return x, np.loadtxt(fname)

class ClsArray(object):
        # Store arrays of cls: self.cls_array[i,j] is zero based array of correlation of field i with j

        def __init__(self, filename=None, cols=None, field_names=['T', 'E', 'B', 'P']):
            self.field_names = field_names
            self.cls_array = np.zeros((len(field_names), len(field_names)), dtype=np.object)
            self.cls_array[:, :] = None
            if filename is not None:
                self.loadFromFile(filename, cols)

        def loadFromFile(self, filename, cols=None):
            if cols is None:
                cols, dat = readWithHeader(filename)
            else:
                dat = np.loadtxt(filename)
            Lix = cols.index('L')
            L = dat[:, Lix]
            self.lmin = L[0]
            self.lmax = L[-1]
            for i, f in enumerate(self.field_names):
                for j, f2 in enumerate(self.field_names[:i + 1]):
                    try:
                        ix = cols.index(f + f2)
                    except:
                        try:
                            ix = cols.index(f2 + f)
                        except:
                            continue
                    cls = np.zeros(self.lmax + 1)
                    cls[self.lmin:self.lmax + 1] = dat[:, ix]
                    self.cls_array[i, j] = cls

        def get(self, indices):
            i, j = indices
            if j > i: i, j = j, i
            return self.cls_array[i, j]

class BinWindows(object):

        def __init__(self, lmin, lmax, nbins):
            self.lmin = lmin
            self.lmax = lmax
            self.nbins = nbins

        def bin(self, TheoryCls, b, cls=None):
            if cls is None: cls = np.zeros(max([x for x in self.cols_out if x is not None]) + 1)
            for i, (pair_in, ix_out) in enumerate(zip(self.cols_in, self.cols_out)):
                cl = TheoryCls.get(pair_in)
                if cl is not None and ix_out is not None:
                    cls[ix_out] += np.dot(self.binning_matrix[b, i, :], cl[self.lmin:self.lmax + 1])
            return cls

        def write(self, froot, stem):
            if not os.path.exists(froot + stem + '_window'): os.mkdir(froot + '_window')
            for b in range(self.nbins):
                with open(froot + stem + '_window/window%u.dat' % (b + 1), 'w') as f:
                    for L in np.arange(self.lmin[b], self.lmax[b] + 1):
                        f.write(("%5u " + "%10e" * len(self.cols_in) + "\n") % (L, self.binning_matrix[b, :, L]))


class DatasetLikelihood(object):

        def __init__(self, fname, field_names=['T', 'E', 'B', 'P']):
            self.field_names = field_names
            self.tot_fields = len(field_names)
            if '.dataset' in fname: self.loadDataset(fname)
            else: raise Exception('lensLike only supports .dataset files')

        def typeIndex(self, field):
            return self.field_names.index(field)

        def clString_to_fieldPair(self, cl):
                if '_' in cl: fields = cl.split('_')
                else:
                    if len(cl) != 2: raise Exception('Cl_order but be CL names, 2 characters or _ separated')
                    fields = [cl[0], cl[1]]
                if len(fields) != 2: raise Exception('Invalid C_l order, must have pairs of field names')
                pair = [self.typeIndex(fields[0]), self.typeIndex(fields[1])]
                if pair[1] > pair[0]: pair.reverse()
                return pair

        def UseString_to_theoryPairCls(self, L):
            pairs = []
            for cl in L:
                pairs.append(self.clString_to_fieldPair(cl))
            return pairs

        def UseString_to_cols(self, L):
            cols = [None] * len(L)
            for i, cl in enumerate(L):
                pair = self.clString_to_fieldPair(cl)
                i1 = self.field_index[pair[0]]
                i2 = self.field_index[pair[1]]

                if i1 < 0 or i2 < 0: continue
                if i2 > i1: i1, i2 = i2, i1
                ix = 0
                for ii in range(self.nfields):
                    for jj in range(ii + 1):
                        if ii == i1 and jj == i2: cols[i] = ix
                        ix += 1
            return cols


        def readBinWindows(self, ini, file_stem):
            bins = BinWindows(self.cl_lmin, self.cl_lmax, self.nbins)
            in_cl = ini.split(file_stem + '_in_order')
            out_cl = ini.split(file_stem + '_out_order')
            bins.cols_in = self.UseString_to_theoryPairCls(in_cl)
            bins.cols_out = self.UseString_to_cols(out_cl)
            norder = len(bins.cols_in)
            if norder != len(bins.cols_out):
                raise Exception('_in_order and _out_order must have same number of entries')

            bins.binning_matrix = np.zeros((self.nbins, norder, self.cl_lmax - self.cl_lmin + 1))
            windows = ini.relativeFileName(file_stem + '_files')
            for b in range(self.nbins):
                window = np.loadtxt(windows % (b + 1))
                Err = False
                for i, L in enumerate(window[:, 0].astype(int)):
                    if self.cl_lmin <= L <= self.cl_lmax:
                        bins.binning_matrix[b, :, L - self.cl_lmin] = window[i, 1:]
                    else:
                        Err = Err or any(window[i, 1:] != 0)
                if Err: print('WARNING: %s %u outside cl_lmin-cl_max range: %s' % (file_stem, b, windows % (b + 1)))
            return bins


        def loadDataset(self, froot):
            if not '.dataset' in froot: froot += '.dataset'
            ini = IniFile(froot)
            self.readIni(ini)

        def readIni(self, ini):
            self.like_approx = ini.string('like_approx', 'gaussian')
            if self.like_approx != 'gaussian': raise Exception('Only gaussian implented in python so far')

            self.fields_use = ini.split('fields_use')
            index_use = [self.typeIndex(f) for f in self.fields_use]
            self.use_field = [i in index_use for i in range(len(self.field_names))]
            self.nfields = sum(self.use_field)

            if ini.hasKey('fields_required'):
                self.fields_required = ini.string('fields_required').split()
            else: self.fields_required = self.fields_use
            index_use = [self.typeIndex(f) for f in self.fields_required]
            self.required_field = [i in index_use for i in range(len(self.field_names))]

            self.binned = ini.bool('binned', True)
            if not self.binned: raise Exception('Currently only support binned')

            self.field_index = np.zeros(self.tot_fields, dtype=int) - 1
            self.fields = np.zeros(self.tot_fields)
            self.field_order = []
            ix = 0
            for i in range(self.tot_fields):
                if self.use_field[i]:
                    self.field_index[i] = ix
                    self.fields[ix] = i
                    self.field_order.append(self.field_names[i])
                    ix += 1
            self.ncl = self.nfields * (self.nfields + 1)

            self.nbins = ini.int('nbins')
            self.bin_min = ini.int('use_min', 1) - 1
            self.bin_max = ini.int('use_max', self.nbins) - 1
            self.nbins_used = self.bin_max - self.bin_min + 1
            self.cl_lmax = ini.int('cl_lmax')
            self.cl_lmin = ini.int('cl_lmin')

            self.phi_lmax = self.cl_lmax
            self.lmin, self.lmax, self.lav, self.bandpowers, self.Ahat = readTextCommentColumns(
                                ini.relativeFileName('cl_hat_file'), ['L_min', 'L_max', 'L_av', 'PP', 'Ahat'])

            self.bins = self.readBinWindows(ini, 'bin_window')

            self.cov = np.loadtxt(ini.relativeFileName('covmat_fiducial'))
            cov = self.cov[self.bin_min:self.bin_max + 1, self.bin_min:self.bin_max + 1]
            self.covinv = np.linalg.inv(cov)

            if 'linear_correction_fiducial_file' in ini.params:
                self.fid_correction = np.loadtxt(ini.relativeFileName('linear_correction_fiducial_file'))[:, 1]
                self.linear_correction = self.readBinWindows(ini, 'linear_correction_bin_window')
            else:
                self.linear_correction = None


        def writeData(self, froot):
            np.savetxt(froot + '_cov.dat', self.cov)
#            self.saveCl(froot + '_fid_cl.dat', self.fid_cl[:, 1:], cols=['TT', 'EE', 'TE', 'PP'])

            with open(froot + '_bandpowers.dat', 'w') as f:
                f.write("#%4s %5s %5s %8s %12s %10s %7s\n" % ('bin', 'L_min', 'L_max', 'L_av', 'PP', 'Error', 'Ahat'))
                for b in range(self.nbins):
                    f.write("%5u %5u %5u %8.2f %12.5e %10.3e %7.3f\n" % (b + 1, self.lmin[b], self.lmax[b], self.lav[b],
                                                         self.bandpowers[b], np.sqrt(self.cov[b, b]), self.Ahat[b]))
            self.bins.write(froot, 'bin')
            if self.linear_correction is not None:
                self.linear_correction.write(froot, 'linear_correction_bin')

            with open(froot + '_lensing_fiducial_correction', 'w') as f:
                f.write("#%4s %12s \n" % ('bin', 'PP'))
                for b in range(self.nbins):
                    f.write("%5u %12.5e\n" % (b + 1, self.fid_correction[b]))

        def plot(self, phicl=None, ls=None):
            lbin = self.lav
            binned_phicl_err = np.zeros(self.nbins)
            for b in range(self.nbins):
                binned_phicl_err[b] = np.sqrt(self.cov[b, b])
            plt.errorbar(lbin, self.bandpowers, yerr=binned_phicl_err  , xerr=[lbin - self.lmin, self.lmax - lbin], fmt='o')

            if phicl is not None:
                if isinstance(phicl, ClsArray): phicl = phicl.get([3, 3])
                if ls is None: ls = np.arange(len(phicl))
                plt.plot(ls, phicl, color='k')
                plt.xlim([2, ls[-1]])

        def chi_squared(self, ClArray):

            binphi = np.zeros(self.nbins_used)
            for b in range(self.bin_min, self.bin_max + 1):
                band = self.bins.bin(ClArray, b)[0]
                if self.linear_correction:
                    band += self.linear_correction.bin(ClArray, b) - self.fid_correction[b]
                binphi[b - self.bin_min] = band

            delta = binphi - self.bandpowers[self.bin_min:self.bin_max + 1]
            return np.dot(delta, np.dot(self.covinv, delta))

def plotAndChisq(dataset, cl_file):
    d = DatasetLikelihood(dataset)
    cls = ClsArray(cl_file)
    d.plot(cls)
    print('Chi-squared: ', d.chi_squared(cls))
    plt.show()


if __name__ == "__main__":
#    plotAndChisq(r'test_data/g60_full_pp.dataset', r'test_data/testout_pp.theory_cl')
#    sys.exit()
    try: import argparse
    except:
        print('use "module load" to load python 2.7')
        sys.exit()
    parser = argparse.ArgumentParser(description="Load .dataset and calculate likelihood")
    parser.add_argument('dataset', help='.dataset filename')
    parser.add_argument('cl_file', help='file of Cls')
    args = parser.parse_args()
    plotAndChisq(args.dataset, args.cl_file)
