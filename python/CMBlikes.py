# Load CosmoMC format .dataset files with lensing likelihood data
# AL July 2014

import numpy as np
import os
import iniFile
from pylab import *

def readTextCommentColumns(fname, cols):
        with open(fname) as f:
            x = f.readline().strip()
            if x[0] != '#': raise Exception('No Comment')
        incols = x[1:].split()
        colnums = [incols.index(col) for col in cols]
        return loadtxt(fname, usecols=colnums, unpack=True)

def readWithHeader(fname):
        with open(fname) as f:
            x = f.readline().strip()
            if x[0] != '#': raise Exception('No Comment')
            x = x[1:].split()
        return x, loadtxt(fname)

class MapCls():
    pass

class ClsArray():
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
                dat = loadtxt(filename)
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

class BinWindows():

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
                        f.write("%5u " + "%10e"*len(self.cols_in) + "\n" % (L, self.binning_matrix[b, :, L]))


class DatasetLikelihood():

        def __init__(self, fname, field_names=['T', 'E', 'B', 'P']):
            self.field_names = field_names
            self.tot_fields = len(field_names)
            if '.dataset' in fname: self.loadDataset(fname)
            else: raise Exception('lensLike only supports .dataset files')

        def typeIndex(self, field):
            return self.field_names.index(field)

        def pairStringToMapIndices(self, S):
            if len(S) == 2:
                if self.has_map_names: raise Exception('CMBlikes: CL names must use MAP1xMAP2 name')
                map1, map2 = S[0], S[1]
            else:
                if not 'x' in S: raise Exception('CMBLikes: invalid spectrum name ' + S)
                map1, map2 = S.split('x')
            return self.map_names.index(map1), self.map_names.index(map2)

        def pairStringToUsedMapIndices(self, used_index, S):
            i1, i2 = self.pairStringToMapIndices(S)
            i1, i2 = used_index[i1], used_index[i2]
            if i2 > i1: i1, i2 = i2, i1
            return i1, i2

        def useString_to_cols(self, S):
            cL_i_j = self.useString_to_Cl_i_j(S, self.map_used_index)

            cols = []
            for i1, i2 in cL_i_j:
                if i1 is None or i2 is None:
                    cols.append(None)
                    continue
                ix = 0
                for ii in range(self.nmaps):
                    for jj in range(ii + 1):
                        if (ii == i1 and jj == i2): cols.append(ix)
                        ix += 1
            return cols

        def useString_to_Cl_i_j(self, S, used_index):
            cl_i_j = []
            if isinstance(S, basestring): S = S.split()
            for P in S:
                cl_i_j.append(self.pairStringToUsedMapIndices(used_index, P))
            return cl_i_j

        def requiredMapPair_to_Theory_i_j(self, i1, i2):
            i = self.map_fields[self.required_order[i1]]
            j = self.map_fields[self.required_order[i2]]
            if j > i: i, j = j, i
            return i, j

        def CL_i_j_name(self, i, j):
            name1, name2 = self.map_order[i], self.map_order[j]
            if self.has_map_names:
                return name1 + 'x' + name2
            else:
                return name1 + name2

        def readBinWindows(self, ini, file_stem):
            bins = BinWindows(self.cl_lmin, self.cl_lmax, self.nbins)
            in_cl = ini.split(file_stem + '_in_order')
            out_cl = ini.split(file_stem + '_out_order')
            bins.cols_in = self.useString_to_Cl_i_j(in_cl, self.map_required_index)
            bins.cols_out = self.useString_to_cols(out_cl)
            norder = len(bins.cols_in)
            if norder != len(bins.cols_out):
                raise Exception('_in_order and _out_order must have same number of entries')

            bins.binning_matrix = np.zeros((self.nbins, norder, self.cl_lmax - self.cl_lmin + 1))
            windows = ini.relativeFileName(file_stem + '_files')
            for b in range(self.nbins):
                window = np.loadtxt(windows % (b + 1))
                Err = False
                for i, L in enumerate(window[:, 0].astype(int)):
                    if (L >= self.cl_lmin and L <= self.cl_lmax):
                        bins.binning_matrix[b, :, L - self.cl_lmin] = window[i, 1:]
                    else:
                        Err = Err or any(window[i, 1:] != 0)
                if Err: print 'WARNING: %s %u outside cl_lmin-cl_max range: %s' % (file_stem, b, windows % (b + 1))

            if ini.hasKey(file_stem + '_fix_cl_file'):
                raise Exception('CMBlikes: python version does not support _fix_cl_file yet')
            return bins


        def loadDataset(self, froot):
            if not '.dataset' in froot: froot += '.dataset'
            ini = iniFile.iniFile(froot)
            self.readIni(ini)

        def readIni(self, ini):
            names = ini.split('map_names', '')
            self.has_map_names = len(names) > 0
            if self.has_map_names:
                self.map_names = names
                self.map_field_names = ini.split('map_fields')
                if len(self.map_field_names) != len(names):
                    raise Exception('CMBLikes: number of map_fields does not match map_names')
                self.map_fields = [self.typeIndex(field) for field in self.map_field_names]
            else:
                self.map_names = self.field_names
                self.map_field_names = self.field_names
                self.map_fields = range(len(self.field_names))

            self.like_approx = ini.string('like_approx', 'gaussian')
            if self.like_approx != 'gaussian': raise Exception('Only gaussian implented in python so far')

            fields_use = ini.split('fields_use', '')
            if len(fields_use):
                use_theory_field = [field in fields_use for field in self.field_names]
            else:
                use_theory_field = [True] * len(self.field_names)

            maps_use = ini.split('maps_use', '')
            if len(maps_use):
                self.use_map = [map in maps_use for map in self.map_names]
            else:
                self.use_map = [use_theory_field[i] for i in self.map_fields]

            if self.has_map_names:
                S = ini.split('maps_required', '')
                if ini.hasKey('fields_required'):
                    raise Exception('CMBLikes: use maps_required not fields_required')
            else:
                S = ini.split('fields_required', '')

            if len(S):
                self.require_map = [m in S for m in self.map_names]
            else:
                self.require_map = self.use_map

            self.required_theory_field = [False] * len(self.field_names)
            for i, r  in enumerate(self.require_map):
                if r:
                    self.required_theory_field[self.map_fields[i]] = True

            self.nmaps = sum(self.use_map)
            self.nmaps_required = sum(self.require_map)

            self.required_order = []
            self.map_required_index = [None] * len(self.map_names)
            self.map_used_index = [None] * len(self.map_names)
            self.map_order = []
            for i, m in enumerate(self.map_names):
                if self.require_map[i]:
                    self.map_required_index[i] = len(self.required_order)
                    self.required_order.append(i)
                if self.use_map[i]:
                    self.map_used_index[i] = len(self.map_order)
                    self.map_order.append(m)

            self.binned = ini.bool('binned', True)
            if not self.binned: raise Exception('Currently only support binned')

            self.ncl = self.nmaps * (self.nmaps + 1)

            self.nbins = ini.int('nbins')
            self.bin_min = ini.int('use_min', 1) - 1
            self.bin_max = ini.int('use_max', self.nbins) - 1
            self.nbins_used = self.bin_max - self.bin_min + 1
            self.cl_lmax = ini.int('cl_lmax')
            self.cl_lmin = ini.int('cl_lmin')

            self.lmin, self.lmax, self.lav, self.bandpowers, self.Ahat = readTextCommentColumns(
                                ini.relativeFileName('cl_hat_file'), ['L_min', 'L_max', 'L_av', 'PP', 'Ahat'])

            self.bins = self.readBinWindows(ini, 'bin_window')

            self.cov = np.loadtxt(ini.relativeFileName('covmat_fiducial'))
            cov = self.cov[self.bin_min:self.bin_max + 1, self.bin_min:self.bin_max + 1]
            self.covinv = inv(cov)

            if 'linear_correction_fiducial' in ini.params:
                self.fid_correction = loadtxt(ini.relativeFileName('linear_correction_fiducial'))[:, 1]
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
            errorbar(lbin, self.bandpowers, yerr=binned_phicl_err  , xerr=[lbin - self.lmin, self.lmax - lbin], fmt='o')
#           if dofid:
#               self.fid_bandpowers = np.dot(self.binning_matrix, self.fid_phi)
#               errorbar(lbin, self.fid_bandpowers , yerr=binned_phicl_err , xerr=[lbin - self.lmin, self.lmax - lbin], fmt='o', alpha=0.3)
#               if phicl is None: phicl = self.fid_phi

            if phicl is not None:
                if isinstance(phicl, ClsArray): phicl = phicl.get([3, 3])
                if ls is None: ls = np.arange(len(phicl))
                plot(ls, phicl, color='k')
                xlim([2, ls[-1]])

        def getTheoryMapCl(self, ClArray):
            used_fields = []
            for name in self.field_names:
                if name in self.map_field_names:used_fields.append(name)
            Cls = ClsArray(field_names=used_fields)
            Cls.lmin = ClArray.lmin
            Cls.lmax = ClArray.lmax
            for i in range(self.nmaps_required):
                for j in range(i + 1):
                    f1, f2 = self.requiredMapPair_to_Theory_i_j(i, j)
                    Cls.cls_array[i, j] = ClArray.cls_array[f1, f2]
            return Cls

        def chi_squared(self, ClArray):

            binphi = np.zeros(self.nbins_used)
            MapCl = self.getTheoryMapCl(ClArray)
            for b in range(self.bin_min, self.bin_max + 1):
                band = self.bins.bin(MapCl, b)[0]
                if self.linear_correction:
                    band += self.linear_correction.bin(MapCl, b) - self.fid_correction[b]
                binphi[b - self.bin_min] = band

            delta = binphi - self.bandpowers[self.bin_min:self.bin_max + 1]
            return np.dot(delta, np.dot(self.covinv, delta))

def plotAndChisq(dataset, cl_file):
    d = DatasetLikelihood(dataset)
    cls = ClsArray(cl_file)
    d.plot(cls)
    print 'Chi-squared: ', d.chi_squared(cls)
    show()


if __name__ == "__main__":
#    plotAndChisq(r'../data/planck_lensing/smica_g30_ftl_full_pp.dataset', r'../data/base_plikHM_TT_lowTEB.minimum.theory_cl')
#    sys.exit()
    try: import argparse
    except:
        print 'use "module load" to load python 2.7'
        sys.exit()
    parser = argparse.ArgumentParser(description="Load .dataset and calculate likelihood")
    parser.add_argument('dataset', help='.dataset filename')
    parser.add_argument('cl_file', help='file of Cls')
    args = parser.parse_args()
    plotAndChisq(args.dataset, args.cl_file)
