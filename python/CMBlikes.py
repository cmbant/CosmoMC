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
            if any([len(f) > 1 for f in field_names]):
                self.sep = 'x'
            else:
                self.sep = ''
            self.cls_array = np.zeros((len(field_names), len(field_names)), dtype=np.object)
            self.cls_array[:, :] = None
            if filename is not None:
                self.loadFromFile(filename, cols)

        def loadFromFile(self, filename, cols=None):
            if cols is None:
                cols, dat = readWithHeader(filename)
            else:
                dat = loadtxt(filename)
                cols = ['L'] + cols
            Lix = 0
            L = dat[:, Lix]
            self.lmin = L[0]
            self.lmax = L[-1]
            for i, f in enumerate(self.field_names):
                for j, f2 in enumerate(self.field_names[:i + 1]):
                    try:
                        ix = cols.index(f + self.sep + f2)
                    except:
                        try:
                            ix = cols.index(f2 + self.sep + f)
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

        def bin_average_l(self):
            avs = np.zeros(self.nbins)
            ls = np.arange(self.lmin, self.lmax + 1)
            for b in range(self.nbins):
                avs[b] = np.dot(self.binning_matrix[b, 0, :], ls) / np.sum(self.binning_matrix[b, 0, :])
            return avs

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
            if i2 > i1:
                return i2, i1
            else:
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

        def mapPair_to_Theory_i_j(self, order, i1, i2):
            i = self.map_fields[order[i1]]
            j = self.map_fields[order[i2]]
            if j > i:
                return j, i
            else:
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
            out_cl = ini.split(file_stem + '_out_order', in_cl)
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
            self.map_required_order = []
            self.map_used_fields = []
            for i, m in enumerate(self.map_names):
                if self.require_map[i]:
                    self.map_required_index[i] = len(self.required_order)
                    self.required_order.append(i)
                    self.map_required_order.append(m)
                if self.use_map[i]:
                    self.map_used_index[i] = len(self.map_order)
                    self.map_order.append(m)
                    self.map_used_fields.append(self.map_fields[i])

            self.binned = ini.bool('binned', True)
            if not self.binned: raise Exception('Currently only support binned')

            self.ncl = (self.nmaps * (self.nmaps + 1)) / 2

            self.nbins = ini.int('nbins')
            self.bin_min = ini.int('use_min', 1) - 1
            self.bin_max = ini.int('use_max', self.nbins) - 1
            self.nbins_used = self.bin_max - self.bin_min + 1
            self.cl_lmax = ini.int('cl_lmax')
            self.cl_lmin = ini.int('cl_lmin')

            self.bandpowers, (self.lmin, self.lmax, self.lav) = self.readClArr(ini, 'cl_hat', info=['L_min', 'L_max', 'L_av'])

#            self.lmin, self.lmax, self.lav, self.bandpowers, self.Ahat = readTextCommentColumns(
#                                ini.relativeFileName('cl_hat_file'), ['L_min', 'L_max', 'L_av', 'PP', 'Ahat'])

            self.cl_hat_includes_noise = ini.bool('cl_hat_includes_noise', False)

            if self.like_approx != 'gaussian' or self.cl_hat_includes_noise:
                self.cl_noise , _ = self.readClArr(ini, 'cl_noise')
                if self.cl_hat_includes_noise:
                    self.bandpowers -= self.cl_noise

            if self.like_approx == 'HL':
                self.cl_fiducial, _ = self.readClArr(ini, 'cl_fiducial')
                self.cl_fiducial_includes_noise = ini.bool('cl_fiducial_includes_noise', False)
                if not self.cl_fiducial_includes_noise:
                    self.cl_fiducial += self.cl_noise
                self.Chat_matrix = []
                self.sqrt_fiducial = []
                for b in range(self.bin_min, self.bin_max + 1):
                    self.Chat_matrix.append(self.elementsToMatrix(self.bandpowers[:, b] + self.cl_noise[:, b]))
                    M = self.elementsToMatrix(self.cl_fiducial[:, b])
                    self.sqrt_fiducial.append(self.matrix_sqrt(M))
            elif self.like_approx == 'exact':
                self.fullsky_exact_fksy = ini.float('fullsky_exact_fksy', 1)

            self.bins = self.readBinWindows(ini, 'bin_window')
            if self.lav is None:
                self.lav = self.bins.bin_average_l()

            self.readCovmat(ini)

            if 'linear_correction_fiducial_file' in ini.params:
                self.fid_correction, _ = self.readClArr(ini, 'linear_correction_fiducial')
                self.linear_correction = self.readBinWindows(ini, 'linear_correction_bin_window')
            else:
                self.linear_correction = None


        def matrix_sqrt(self, M):
            diag, U = np.linalg.eigh(M)
            diag = np.sqrt(diag)
            return (diag * U).dot(U.T)

        def readCovmat(self, ini):
            cov_cl = ini.split('covmat_cl')
            cl_in_index = self.useString_to_cols(cov_cl)
            num_in = len(cl_in_index)
            self.ncl_used = sum([True for x in cl_in_index if x is not None])
            ix = -1
            self.cl_use_index = np.ndarray(self.ncl_used, dtype=int)
            cov_cl_used = np.ndarray(self.ncl_used, dtype=int)
            for i, index in enumerate(cl_in_index):
                if index is not None:
                    ix += 1
                    self.cl_use_index[ix] = index
                    cov_cl_used[ix] = i

            self.cov = np.loadtxt(ini.relativeFileName('covmat_fiducial'))
            cov = np.zeros((self.nbins_used * self.ncl_used, self.nbins_used * self.ncl_used))
            for binx in range(self.bin_min, self.bin_max + 1):
                for biny in range(self.bin_min, self.bin_max + 1):
                    cov[(binx - self.bin_min) * self.ncl_used:(binx - self.bin_min + 1) * self.ncl_used,
                        (biny - self.bin_min) * self.ncl_used:(biny - self.bin_min + 1) * self.ncl_used ] = \
                        self.cov[np.ix_(binx * num_in + cov_cl_used, biny * num_in + cov_cl_used)]
#            cov = self.cov[self.bin_min:self.bin_max + 1, self.bin_min:self.bin_max + 1]
            self.covinv = inv(cov)


#         def writeData(self, froot):
#             np.savetxt(froot + '_cov.dat', self.cov)
# #            self.saveCl(froot + '_fid_cl.dat', self.fid_cl[:, 1:], cols=['TT', 'EE', 'TE', 'PP'])
#
#             with open(froot + '_bandpowers.dat', 'w') as f:
#                 f.write("#%4s %5s %5s %8s %12s %10s %7s\n" % ('bin', 'L_min', 'L_max', 'L_av', 'PP', 'Error', 'Ahat'))
#                 for b in range(self.nbins):
#                     f.write("%5u %5u %5u %8.2f %12.5e %10.3e %7.3f\n" % (b + 1, self.lmin[b], self.lmax[b], self.lav[b],
#                                                          self.bandpowers[b], np.sqrt(self.cov[b, b]), self.Ahat[b]))
#             self.bins.write(froot, 'bin')
#             if self.linear_correction is not None:
#                 self.linear_correction.write(froot, 'linear_correction_bin')
#
#             with open(froot + '_lensing_fiducial_correction', 'w') as f:
#                 f.write("#%4s %12s \n" % ('bin', 'PP'))
#                 for b in range(self.nbins):
#                     f.write("%5u %12.5e\n" % (b + 1, self.fid_correction[b]))

        def Cl_i_j_name(self, name1, name2):
            if self.has_map_names:
                return name1 + 'x' + name2
            else:
                return name1 + name2

        def readClArr(self, ini, basename, info=[]):
            filename = ini.relativeFileName(basename + '_file')
            order = ini.split(basename + '_order', '')
            if not order:
                order, dat = readWithHeader(filename)
            else:
                dat = loadtxt(filename)
                order = ['bin'] + order
            cols = np.zeros((self.ncl, dat.shape[0]))

            if dat[0, 0] != 1: raise Exception('Expect file of bin numbers and values')
            if dat.shape[0] < self.bin_max: raise Exception('file does not go up to maximum bin:' + filename)
            ext_info = []
            for ext in info:
                try:
                    cix = order.index(ext)
                    ext_info.append(dat[:, cix])
                except:
                    ext_info.append(None)
                    continue
            ix = -1
            for i, f in enumerate(self.map_order):
                for f2 in self.map_order[:i + 1]:
                    ix += 1
                    try:
                        cix = order.index(self.Cl_i_j_name(f, f2))
                    except:
                        try:
                            cix = order.index(self.Cl_i_j_name(f2 , f))
                        except:
                            continue
                    cols[ix, :] = dat[:, cix]
            return cols, ext_info

        def elementsToMatrix(self, X, outM=None):
            M = outM or np.ndarray((self.nmaps, self.nmaps), dtype=X.dtype)
            off = 0
            for i in range(self.nmaps):
                M[i, 0:i + 1] = X[off:off + i + 1]
                M[0:i + 1, i] = X[off:off + i + 1]
                off += i + 1
            return M

        def matrixToElements(self, M, outX=None):
            X = outX or np.ndarray(self.ncl, dtype=M.dtype)
            off = 0
            for i in range(self.nmaps):
                X[off:off + i + 1] = M[i, 0:i + 1]
                off += i + 1
            return X

        def required_i_j_for_used_fields(self, cl, i, j):
            f1, f2 = self.map_used_fields[i], self.map_used_fields[j]
            for i1 in range(self.nmaps_required):
                for j1 in range(i1 + 1):
                    ff1, ff2 = cl.theory_i_j[i1, j1]
                    if ff1 == f1 and ff2 == f2:
                        return i1, j1


        def plotBands(self, cl=None, lmax=None):
            lbin = self.lav
            binned_phicl_err = np.zeros(self.nbins)
            if self.nmaps > 1: f, plotax = plt.subplots(self.nmaps, self.nmaps)
            ix = -1
            for i in range(self.nmaps):
                for j in range(self.nmaps):
                    if j > i:
                        plotax[i, j].set_visible(False)
                        continue
                    ix += 1
                    if not ix in self.cl_use_index: continue
                    if self.nmaps > 1:
                        ax = plotax[i, j]
                        ax.set_title(self.CL_i_j_name(i, j))
                    else: ax = gca()
                    this_ix = list(self.cl_use_index).index(ix)
                    for b in range(self.nbins):
                        binned_phicl_err[b] = np.sqrt(self.cov[b * self.ncl_used + this_ix, b * self.ncl_used + this_ix])

                    if self.lmin is not None:
                        ax.errorbar(lbin, self.bandpowers[ix, :], yerr=binned_phicl_err  , xerr=[lbin - self.lmin, self.lmax - lbin], fmt='o')
                    else:
                        ax.errorbar(lbin, self.bandpowers[ix, :], yerr=binned_phicl_err , fmt='o')

                    if cl is not None:
                        if self.nmaps_required > self.nmaps:
                            i1, j1 = self.required_i_j_for_used_fields(cl, i, j)
                        else:
                            i1, j1 = i, j
                        lmax = lmax or min(cl.lmax, self.cl_lmax)
                        ls = np.arange(2, lmax + 1)
                        acl = cl.get((i1, j1))
                        ax.plot(ls, acl[2:lmax + 1], color='k')
                        ax.set_xlim([2, lmax])

#           if dofid:
#               self.fid_bandpowers = np.dot(self.binning_matrix, self.fid_phi)
#               errorbar(lbin, self.fid_bandpowers , yerr=binned_phicl_err , xerr=[lbin - self.lmin, self.lmax - lbin], fmt='o', alpha=0.3)
#               if phicl is None: phicl = self.fid_phi


        def getTheoryMapCl(self, ClArray):
            Cls = ClsArray(field_names=self.map_required_order)
            Cls.lmin = ClArray.lmin
            Cls.lmax = ClArray.lmax
            Cls.theory_i_j = np.ndarray((self.nmaps_required, self.nmaps_required), dtype=object)
            for i in range(self.nmaps_required):
                for j in range(i + 1):
                    f1, f2 = self.mapPair_to_Theory_i_j(self.required_order, i, j)
                    Cls.theory_i_j[i, j] = f1, f2
                    if ClArray.cls_array[f1, f2] is not None:
                        Cls.cls_array[i, j] = np.copy(ClArray.cls_array[f1, f2])
            return Cls

        def getBinnedMapCls(self, MapCl, b):
            band = self.bins.bin(MapCl, b)
            if self.linear_correction:
                band += self.linear_correction.bin(MapCl, b) - self.fid_correction[:, b]
            return band

        def HL_transform(self, vecp, Chat, CfHalf):
            C = self.elementsToMatrix(vecp)
            diag, U = np.linalg.eigh(C)  # c = (U * diag).dot(U.T)
            Rot = U.T.dot(Chat).dot(U)
            roots = np.sqrt(diag)
            for i in range(self.nmaps):
                Rot[i, :] /= roots
                Rot[:, i] /= roots
            Rot = U.dot(Rot).dot(U.T)
            diag, Rot = np.linalg.eigh(Rot)
            tmp = diag - np.log(diag) - 1
            tmp[tmp < 0] = 0
            diag = np.sign(diag - 1) * (np.sqrt(2 * tmp))

            U = CfHalf.dot(Rot)
            C = np.copy(U)
            for i in range(self.nmaps):
                C[:, i] *= diag[i]
            C = C.dot(U.T)
            return self.matrixToElements(C)


        def loadForegroundFile(self, filename):
            return ClsArray(filename, field_names=self.map_order)

        def addForegrounds(self, MapCl, ForegroundCl):
            lmax = min(MapCl.lmax, ForegroundCl.lmax)
            for i in range(self.nmaps_required):
                for j in range(i + 1):
                    if MapCl.cls_array[i, j] is None: MapCl.cls_array[i, j] = ForegroundCl.cls_array[i, j]
                    elif ForegroundCl.cls_array[i, j] is not None:
                        MapCl.cls_array[i, j][:lmax] += ForegroundCl.cls_array[i, j][:lmax]

        def chi_squared(self, MapClArray):
            bigX = np.empty((self.bin_max - self.bin_min + 1) * self.ncl_used)
            for b in range(self.bin_min, self.bin_max + 1):
                boff = b - self.bin_min
                vecp = self.getBinnedMapCls(MapClArray, b)
                if self.like_approx == 'HL':
                    vecp = vecp + self.cl_noise[:, b]
                    vecp = self.HL_transform(vecp, self.Chat_matrix[boff], self.sqrt_fiducial[boff])
                else:
                    vecp -= self.bandpowers[:, b]
                bigX[boff * self.ncl_used:(boff + 1) * self.ncl_used] = vecp[self.cl_use_index]

            return np.dot(bigX, np.dot(self.covinv, bigX))


def plotAndChisq(dataset, cl_file, foreground_file=None, lmax=None):
    d = DatasetLikelihood(dataset)
    cls = ClsArray(cl_file)
    MapCl = d.getTheoryMapCl(cls)
    if foreground_file:
        foregrounds = d.loadForegroundFile(foreground_file)
        d.addForegrounds(MapCl, foregrounds)
    d.plotBands(MapCl, lmax)
    print 'Chi-squared: ', d.chi_squared(MapCl)
    show()


if __name__ == "__main__":
#    plotAndChisq(r'../data/planck_lensing/smica_g30_ftl_full_pp.dataset', r'../data/base_plikHM_TT_lowTEB.minimum.theory_cl')
#   plotAndChisq(r'../data/bicep/BK_Planck.dataset', r'z://testout_pp.theory_cl', r'z://testout_pp.BKPLANCK_foregrounds', 200)
#   sys.exit()
    try: import argparse
    except:
        print 'use "module load" to load python 2.7'
        sys.exit()
    parser = argparse.ArgumentParser(description="Load .dataset and calculate likelihood")
    parser.add_argument('dataset', help='.dataset filename')
    parser.add_argument('cl_file', help='file of Cls')
    parser.add_argument('--foreground_file', help='file of foreground Cls')
    parser.add_argument('--lmax')

    args = parser.parse_args()
    plotAndChisq(args.dataset, args.cl_file, args.foreground_file, args.lmax)
