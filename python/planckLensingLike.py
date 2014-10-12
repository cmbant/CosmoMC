# Load Duncan lenslike .dat files, convert to cosmomc format
# and calculate likelihood using first order perturbative renormalization corrections
# AL May 2014

from pylab import *
import numpy as np
import string
import os
import iniFile
import CMBlikes

def parseStr(x):
    if x.isdigit(): return int(x)
    if x.isalpha() or x.isalnum(): return x
    return float(x)

def readParsedLine(f, tp=None):
    if tp is None:
        return [parseStr(x) for x in f.readline().strip().split()]
    else:
        return [tp(x) for x in f.readline().strip().split()]


class quad_estimator: pass

def normPhi(L):
        res = (L * (L + 1.0)) ** 2 / np.pi / 2
        if L[0] == 0: res[0] = 1
        return res

def normCMB(L):
        res = (L * (L + 1.0)) / np.pi / 2
        if L[0] == 0: res[0] = 1
        return res

class lensLike(CMBlikes.DatasetLikelihood):

        def __init__(self, fname, bin_compression=False):
            self.bin_compression = bin_compression
            if '.dataset' in fname: self.loadDataset(fname)
            else: self.convertDuncan(fname)

        def convertDuncan(self, fname):
            with open(fname) as f:
                x = '#'
                while x[0] == '#':
                    x = f.readline().strip()
                self.nbins = int(x)
                self.cmb_lmax, self.phi_lmax = readParsedLine(f)
                n1qbins, n1pbins = readParsedLine(f)
                self.num_estimators, self.num_estimator_crosses = readParsedLine(f)
                self.s4hat, self.s4std = readParsedLine(f)
                print '%u nbins, %u estimators' % (self.nbins, self.num_estimators)
                self.lmin = np.zeros(self.nbins)
                self.lmax = np.zeros(self.nbins)
                self.Ahat = np.zeros(self.nbins)
                for ix in range(self.nbins):
                    _, self.lmin[ix], self.lmax[ix], self.Ahat[ix] = readParsedLine(f)
                self.lmid = (self.lmax + self.lmin) / 2
                print 'L midpoints:', self.lmid

                self.cov = np.zeros((self.nbins, self.nbins))
                for ix in range(self.nbins):
                    self.cov[:, ix] = readParsedLine(f)
                self.covinv = np.zeros((self.nbins, self.nbins))
                for ix in range(self.nbins):
                    self.covinv[:, ix] = readParsedLine(f)
                self.fid_cl = np.zeros((self.cmb_lmax + 1, 5))
                for ix in range(self.cmb_lmax + 1):
                    self.fid_cl[ix, :4] = readParsedLine(f)
                for ix in range(1, 4):
                    self.fid_cl[:, ix] *= normCMB(self.fid_cl[:, 0])

                self.fid_phi = np.zeros((self.phi_lmax + 1, 7))
                for ix in range(self.phi_lmax + 1):
                    self.fid_phi[ix, :] = readParsedLine(f)
                self.fid_norm = self.fid_phi[:, 6]

                self.fid_cl = self.fid_cl[0:self.phi_lmax + 1, :]
                self.fid_cl[:, 4] = self.fid_phi[:, 1] * normPhi(self.fid_cl[:, 0])

                self.binning_matrix = np.zeros((self.nbins, self.phi_lmax + 1))
                self.fid_bandpowers = np.zeros(self.nbins)
                self.lav = np.zeros(self.nbins)
                for b in range(self.nbins):
                    phicl = self.fid_phi[self.lmin[b]:self.lmax[b] + 1, 1]
                    tmp = phicl * self.fid_phi[self.lmin[b]:self.lmax[b] + 1, 5]
                    # switch to operate on L(L+1)CL/2pi (propto kappa power), rather than phi
                    # choice of how to bormalize bandpowers fairly abtirary
                    ls = np.arange(self.lmin[b], self.lmax[b] + 1)
                    norm = (ls * (ls + 1.0)) ** 2 / 2 / np.pi
                    tmp = tmp / norm
#                    bnorm = sum(tmp * phicl)
                    bnorm = sum(tmp)
                    tmp /= bnorm

                    self.binning_matrix[b, self.lmin[b]:self.lmax[b] + 1] = tmp

                    # self.lav[b] = sum(tmp * ls * phicl)
                    # a = self.lav[b]
                    self.lav[b] = sum(tmp * ls)
                    print self.lmid[b], self.lav[b]

                    self.fid_bandpowers[b] = sum(tmp * norm * phicl)
#                    self.bandpowers[b] = self.Ahat[b] * self.fid_phi[self.lmid[b], 1] * sum(tmp * phicl)
                self.bandpowers = self.fid_bandpowers * self.Ahat
                for b in range(self.nbins):
                    for b2 in range(self.nbins):
                        self.cov[b, b2] *= self.fid_bandpowers[b] * self.fid_bandpowers[b2]
                        self.covinv[b, b2] /= self.fid_bandpowers[b] * self.fid_bandpowers[b2]

                self.N1_q_L = np.zeros(n1qbins, np.int)
                for i in range(n1qbins):
                    _, self.N1_q_L[i] = readParsedLine(f, tp=int)
                self.N1_p_L = np.zeros(n1pbins, np.int)
                for i in range(n1pbins):
                    _, self.N1_p_L[i] = readParsedLine(f, tp=int)
                self.N1_matrix = np.zeros((n1qbins, n1pbins))
                norm = normPhi(self.N1_p_L)
                normL = normPhi(self.N1_q_L)
                print self.N1_p_L
                for i, L in enumerate(self.N1_q_L):
                    self.N1_matrix[i, :] = normL[i] * np.array(readParsedLine(f)) * self.fid_phi[L, 6] / norm
                print 'N1 shape:', self.N1_matrix.shape
                ls = self.fid_phi[:, 0]
                self.fid_N1 = self.fid_phi[:, 6] * self.fid_phi[:, 4] * normPhi(ls)

                self.fid_phi = self.fid_cl[:, 4]

                # read description of quadratic estimators
                self.estimators = []
                for i in range(self.num_estimators):
                    e = quad_estimator()
                    self.estimators.append(e)
                    e.num_separable_terms, = readParsedLine(f)  # ntrm in c code
                    e.lmax, = readParsedLine(f)
                    e.s12 = []

                    # spins and weights, as documented in the table of the 2014 paper
                    # read spins of the estimators
                    for j in range(e.num_separable_terms):
                        e.s12.append(readParsedLine(f))  # triplet

                    # read weights as function of L for each term
                    e.w12 = []
                    for j in range(e.num_separable_terms):
                        w12L = np.zeros((e.lmax + 1, 3))
                        e.w12.append(w12L)
                        for L in range(e.lmax + 1):
                            tmp = readParsedLine(f)
                            if tmp[0] != L: raise Exception('Error reading weights')
                            w12L[L, :] = tmp[1:]
#                            if i == 0 and j == 0: print L, -2 * tmp[2]

                # read fields
                for i in range(self.num_estimators):
                    tmp, a, b = readParsedLine(f)
                    if tmp != i: raise Exception('Error reading fields')
                    self.estimators[i].fields = [a, b]
                self.crosses = []
                for i in range(self.num_estimator_crosses):
                    tmp, q12, q34 = readParsedLine(f)
                    if tmp != i: raise Exception('Error reading q12, q34')
                    self.crosses.append([q12, q34])
                self.long_hash = f.readline().strip()


        def saveCl(self, fname, arr, cols=None, norm=1):
            with open(fname, 'w') as f:
                if cols is not None:
                    f.write('#%4s ' % ('L') + " ".join(["%12s " % (cl) for cl in cols]) + "\n")
                for L in range(arr.shape[0]):
                    if norm == 'lensing':
                        sc = (L * (L + 1.0)) ** 2 / np.pi / 2
                    if norm == 'cl':
                        sc = (L * (L + 1.0)) / np.pi / 2
                    else:
                        sc = norm
                    if len(arr.shape) == 1:
                        f.write('%5u %12e\n' % (L, arr[L] * sc))
                    else:
                        f.write('%5u ' % (L) + " ".join(["%12e " % (cl) for cl in (arr[L, :] * sc)]) + "\n")

        def loadIndexedMatrix(self, fname):
            M = loadtxt(fname)
            return M[1:, 0].astype(int), M[0, 1:].astype(int), M[1:, 1:]


        def getRenormDataFromMatrix(self, fnameAL, fnameN1):
            # L by l matrix mapping TT to normalization
            # this file starts at L=2, and first column is L
            # bin_renorm_matrix starts at zero as others
            self.renorm_L, self.renorm_LT, self.renorm_matrix = self.loadIndexedMatrix(fnameAL)

            self.renorm_fid = np.dot(self.renorm_matrix, self.fid_cl[self.renorm_LT, 1])
            self.renorm_N1_L, self.renorm_N1_LT, self.renorm_N1_matrix = self.loadIndexedMatrix(fnameN1)

            self.renorm_N1_fid = np.dot(self.renorm_N1_matrix, self.fid_cl[self.renorm_N1_LT, 1])

            # bin_renorm_matrix gets correction to A_L (so square it for effect on CPhi)
#            print self.renorm_L
#            np.multiply(self.fid_phi[2:self.phi_lmax + 1] , self.renorm_matrix[:self.phi_lmax - 1, :])
            if self.bin_compression:
                self.bin_renorm_matrix = np.dot(self.binning_matrix[:, 2:] * self.fid_phi[2:self.phi_lmax + 1], self.renorm_matrix[:self.phi_lmax - 1, :])
                self.bin_renorm_N1_matrix = np.dot(self.binning_matrix[:, 2:], self.renorm_N1_matrix[:self.phi_lmax - 1, :])
                self.bin_N1_matrix = np.dot(self.binning_matrix[:, 2:], self.N1_matrix[:self.phi_lmax - 1, :])
                self.T_window = 2 * self.bin_renorm_N1_matrix + 2 * self.bin_renorm_matrix
                self.fid_correction = np.dot(self.T_window , self.fid_cl[self.renorm_LT, 1])
            #    tmp = np.zeros(self.N1_p_L[-1] + 1)
            #    tmp[:self.fid_phi.shape[0]] = self.fid_phi
                self.fid_correction += np.dot(self.bin_N1_matrix, self.fid_phi[self.N1_p_L])

        def dumpData(self, froot, extra=False):
            # these not needed for likelihood
            if hasattr(self, "estimators") and extra:
                self.saveCl(froot + '_TT_Filter.dat', -2 * self.estimators[0].w12[0][:, 1], cols=['F_l'])

            # main files
            np.savetxt(froot + '_cov.dat', self.cov)
#            self.saveCl(froot + '_fid_cl.dat', self.fid_cl[:, 1:], cols=['TT', 'EE', 'TE', 'PP'])
#           self.saveCl(froot + '_fid_N1.dat', self.fid_N1, cols=['PP'])

            with open(froot + '_bandpowers.dat', 'w') as f:
                f.write("#%4s %5s %5s %8s %12s %10s %7s\n" % ('bin', 'L_min', 'L_max', 'L_av', 'PP', 'Error', 'Ahat'))
                for b in range(self.nbins):
                    f.write("%5u %5u %5u %8.2f %12.5e %10.3e %7.3f\n" % (b + 1, self.lmin[b], self.lmax[b], self.lav[b],
                                                         self.bandpowers[b], np.sqrt(self.cov[b, b]), self.Ahat[b]))
            if not os.path.exists(froot + '_window'): os.mkdir(froot + '_window')
            for b in range(self.nbins):
                with open(froot + '_window/window%u.dat' % (b + 1), 'w') as f:
                    for L in np.arange(self.lmin[b], self.lmax[b] + 1):
                        f.write("%5u %10e\n" % (L, self.binning_matrix[b, L]))

            return

            with open(froot + '_fid_N1_dphi.dat', 'w') as f:
                f.write("%5u " % (0) + "".join("%15u " % (L) for L in self.N1_p_L) + "\n")
                for i, L in enumerate(self.N1_q_L):
                    f.write("%5u " % L + " ".join("%15.8e" % (L2) for L2 in self.N1_matrix[i, :]) + "\n")

            if self.bin_compression:
                Ls = [L for L in self.renorm_LT]
                for L in self.N1_p_L:
                    if not L in Ls: Ls.append(L)
                if not os.path.exists(froot + '_lens_delta_window'): os.mkdir(froot + '_lens_delta_window')
                for b in range(self.nbins):
                    with open(froot + '_lens_delta_window/window%u.dat' % (b + 1), 'w') as f:
                        for L in Ls:
                            if L in self.renorm_LT: T = self.T_window[b, np.where(self.renorm_LT == L)]
                            else: T = 0
                            if L in self.N1_p_L: phi = self.bin_N1_matrix[b, np.where(self.N1_p_L == L)]
                            else: phi = 0
                            f.write("%5u %10e %10e\n" % (L, T, phi))
                with open(froot + '_lensing_fiducial_correction.dat', 'w') as f:
                    f.write("#%4s %12s \n" % ('bin', 'PP'))
                    for b in range(self.nbins):
                        f.write("%5u %12.5e\n" % (b + 1, self.fid_correction[b]))


        def testLikes(self):
            print 'chi2 fid=', lens.chi_squared(self.fid_phi, None)
            ls = self.fid_cl[:2049, 0]
            cl_fid = self.fid_cl[:2049, 1]
            cl2 = cl_fid * (ls / 1000.0) ** 0.03 * 1.02
            phi = self.fid_phi * 1.1 * (ls / 1000.0) ** (0.02)
            print 'chi2 =', self.chi_squared(phi, cl2)


def ConvertDuncan(fname, dir, root, var):
        lens = lensLike(fname + '.dat', bin_compression=True)
        outroot = root
        if var: outroot += '_' + var
        diroutroot = os.path.join(r'z://lens/', outroot)
        lens.dumpData(diroutroot)
        for b in range(lens.nbins):
            win = loadtxt(fname + '_lens_delta_window/window%u.dat' % (b + 1))
            if not os.path.exists(diroutroot + '_lens_delta_window'): os.mkdir(diroutroot + '_lens_delta_window')
            with open(diroutroot + '_lens_delta_window/window%u.dat' % (b + 1), 'w') as f:
                for i in range(len(win[:, 0])):
                    if lens.num_estimators == 1:
                        f.write("%5u %10e %10e\n" % (int(win[i, 0]), win[i, 1], win[i, 4]))
                    else:
                        f.write("%5u %10e %10e %10e %10e\n" % (int(win[i, 0]), win[i, 1], win[i, 2], win[i, 3], win[i, 4]))
        corr = loadtxt(fname + '_lensing_fiducial_correction.dat')
        with open(diroutroot + '_lensing_fiducial_correction.dat', 'w') as f:
            f.write("#%4s %14s \n" % ('bin', 'PP'))
            for b in range(lens.nbins):
                f.write("%5u %14.7e\n" % (b + 1, corr[b, 1]))


        with open(diroutroot + '.dataset', "w") as f:
            f.write("like_approx = gaussian\n")
            f.write("fields_use = P\n")
            if lens.num_estimators == 1:
                f.write("fields_required =  T\n")
                order = 'TT PP'
                out_order = 'PP PP'
            else:
                f.write("fields_required =  T E P\n")
                order = 'TT EE TE PP'
                out_order = ' PP PP PP PP'
            f.write("binned = T\n")
            f.write("nbins= %u \n" % (lens.nbins))
            f.write("use_min= 1 \n")
            f.write("use_max= %u \n" % (lens.nbins))
            f.write("cl_lmin = 2\n")
            f.write("cl_lmax = %s\n" % (lens.phi_lmax))
            f.write("""
cl_hat_file = %s_bandpowers.dat

bin_window_files = %s_window/window%%u.dat
bin_window_in_order = PP
bin_window_out_order = PP

covmat_cl = PP
covmat_fiducial = %s_cov.dat

linear_correction_fiducial = %s_lensing_fiducial_correction.dat
linear_correction_bin_window_files = %s_lens_delta_window/window%%u.dat
linear_correction_bin_window_in_order = %s
linear_correction_bin_window_out_order = %s
            """ % (outroot, outroot, outroot, outroot, outroot, order, out_order))

            with open(diroutroot + '_lensonly.dataset', 'w') as f:
                    f.write(
"""DEFAULT(%s.dataset)
INCLUDE(fix_renormalization.ini)""" % (outroot))


if __name__ == "__main__":
    dir = r'C:\Work\F90\LensingBiases\plenslike_dx11'
    for root in ['smica_g30_ftl_full_pp', 'smica_g30_ftl_full_pttptt']:
        for var, newvar in [(r'aggressive', 'aggressive'), (r'conservative', '')]:
            Duncan = True
            if Duncan:
                ConvertDuncan(os.path.join(dir , var , root), dir, root, newvar)

            else:
                f = dir + os.sep + root + '.dataset'
                lens = CMBlikes.DatasetLikelihood(f)
        #        lens.dumpData(r'C:\Work\F90\LensingBiases\like' + os.sep + root)

                lens.plot()
                show()

