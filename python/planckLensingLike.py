# Load Duncan lenslike .dat files, convert to cosmomc format
# and calculate likelihood using first order perturbative renormalization corrections
# AL May 2014

from pylab import *
import numpy as np
import string
import os
import iniFile
from scipy import interpolate


def parseStr(x):
    if x.isdigit(): return int(x)
    if x.isalpha() or x.isalnum(): return x
    return float(x)

def readParsedLine(f, tp=None):
    if tp is None:
        return [parseStr(x) for x in f.readline().strip().split()]
    else:
        return [tp(x) for x in f.readline().strip().split()]

def readTextCommentColumns(fname, cols):
        with open(fname) as f:
            x = f.readline().strip()
            if x[0] != '#': raise Exception('No Comment')
        incols = x[1:].split()
        colnums = [incols.index(col) for col in cols]
        return loadtxt(fname, usecols=colnums, unpack=True)

class quad_estimator: pass

def normPhi(L):
        res = (L * (L + 1.0)) ** 2 / np.pi / 2
        if L[0] == 0: res[0] = 1
        return res

def normCMB(L):
        res = (L * (L + 1.0)) / np.pi / 2
        if L[0] == 0: res[0] = 1
        return res

class lensLike():

        def __init__(self, fname):
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


        def makeN1(self, phicl, clTT=None):
            # first get N1 for fiducial ClTT, then rescale result by change in fiducial N1 for actual clTT at fiducial CPhi
            N1samps = np.dot(self.N1_matrix, phicl[self.N1_p_L])
            f = interpolate.interp1d(self.N1_q_L, N1samps, bounds_error=False)
            N1 = f(np.arange(self.phi_lmax + 1))
            if clTT is not None:
#                N1_renorm = (np.dot(self.renorm_N1_matrix, clTT[self.renorm_LT]) / self.renorm_N1_fid - 1)
                N1_renorm = np.dot(self.renorm_N1_matrix, clTT[self.renorm_LT] - self.fid_cl[self.renorm_LT, 1])
                f = interpolate.interp1d(self.renorm_L, N1_renorm, bounds_error=False)
#               print  (2 * f(np.arange(2, self.phi_lmax + 1)))[::20]

#                print N1[2::20] * (2 * f(np.arange(2, self.phi_lmax + 1)))[::20]
#                N1[2:] = N1[2:] * (1 + 2 * f(np.arange(2, self.phi_lmax + 1)))
                N1[2:] += 2 * f(np.arange(2, self.phi_lmax + 1))
                # print 2 * f(np.arange(2, self.phi_lmax + 1, 50))
            return N1

        def makeRenormCPhi(self, phicl, clTT=None):
            if clTT is None: return phicl
            if False:
                renorm = (np.dot(self.renorm_matrix, clTT[self.renorm_LT]) / self.renorm_fid) ** 2
                # print renorm
                f = interpolate.interp1d(self.renorm_L, renorm, bounds_error=False)
                phi = phicl
                phi[2:] = phi[2:] * f(np.arange(2, self.phi_lmax + 1))
            else:
                renorm = np.dot(self.renorm_matrix, clTT[self.renorm_LT] - self.fid_cl[self.renorm_LT, 1])
                # print 1 + 2 * renorm
                f = interpolate.interp1d(self.renorm_L, renorm, bounds_error=False)
                phi = phicl
                phi[2:] = phi[2:] * (1 + 2 * f(np.arange(2, self.phi_lmax + 1)))
            return phi

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
#            tmp = self.renorm_matrix[:self.phi_lmax - 1, 1:]
#            for i in range(tmp.shape[0]):
#                tmp[i, :] = tmp[i, :] / fid_norm[i]
#            for b in range(self.nbins):
#                self.bin_renorm_matrix[b, 2:] = np.dot(self.binning_matrix[b, 2:], tmp)

#            tmp = self.renorm_N1_matrix[:self.phi_lmax - 1, 1:]
#            for i in range(tmp.shape[0]):
#                tmp[i, :] = tmp[i, :] / fid_norm[i]
#            self.bin_N1_renorm_matrix = np.zeros((self.nbins, lmax_mat + 1))
#            for b in range(self.nbins):
#                self.bin_N1_renorm_matrix[b, 2:] = np.dot(self.binning_matrix[b, 2:], tmp)


        def dumpData(self, froot):
            # these not needed for likelihood
            self.saveCl(froot + '_TT_Filter.dat', -2 * self.estimators[0].w12[0][:, 1], cols=['F_l'])
            self.saveCl(froot + '_fid_N1.dat', self.fid_N1, cols=['PP'])

            # main files
            np.savetxt(froot + '_cov.dat', self.cov)
            self.saveCl(froot + '_fid_cl.dat', self.fid_cl[:, 1:], cols=['TT', 'EE', 'TE', 'PP'])

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

            with open(froot + '_fid_N1_dphi.dat', 'w') as f:
                f.write("%5u " % (0) + "".join("%15u " % (L) for L in self.N1_p_L) + "\n")
                for i, L in enumerate(self.N1_q_L):
                    f.write("%5u " % L + " ".join("%15.8e" % (L2) for L2 in self.N1_matrix[i, :]) + "\n")


        def loadDataset(self, froot):
            if not '.dataset' in froot: froot += '.dataset'
            ini = iniFile.iniFile(froot)
            self.nbins = ini.int('nbins')
            self.bin_min = ini.int('use_min', 1) - 1
            self.bin_max = ini.int('use_max', self.nbins) - 1
            self.cmb_lmax = ini.int('cl_lmax')
            self.phi_lmax = self.cmb_lmax
            self.lmin, self.lmax, self.lav, self.bandpowers, self.Ahat = readTextCommentColumns(
                                ini.relativeFileName('cl_hat_file'), ['L_min', 'L_max', 'L_av', 'PP', 'Ahat'])

            self.binning_matrix = np.zeros((self.nbins, self.phi_lmax + 1))
            windows = ini.relativeFileName('bin_window_files')
            for b in range(self.nbins):
                bins = loadtxt(windows % (b + 1))
                self.binning_matrix[b, bins[:, 0].astype(int)] = bins[:, 1]
            self.cov = loadtxt(ini.relativeFileName('covmat_fiducial'))
            self.covinv = inv(self.cov)

            self.fid_cl = loadtxt(ini.relativeFileName('lensing_fiducial_cl'))
            self.fid_phi = self.fid_cl[:, 4]

            self.N1_q_L, self.N1_p_L, self.N1_matrix = self.loadIndexedMatrix(ini.relativeFileName('lensing_N1_matrix_fiducial_dphi'))

            self.getRenormDataFromMatrix(ini.relativeFileName('lensing_renorm_matrix') % ('TT'),
                                          ini.relativeFileName('lensing_renorm_N1_matrix') % ('TT'))

        def plot(self, phicl=None, dofid=True, doN1=True):
            ls = np.arange(self.phi_lmax + 1)
            if phicl is None: phicl = self.fid_phi
            lbin = self.lav
            binned_phicl_err = np.zeros(self.nbins)
            for b in range(self.nbins):
                binned_phicl_err[b] = np.sqrt(self.cov[b, b])
            errorbar(lbin, self.bandpowers, yerr=binned_phicl_err  , xerr=[lbin - self.lmin, self.lmax - lbin], fmt='o')
            if dofid:
                self.fid_bandpowers = np.dot(self.binning_matrix, self.fid_phi)
                errorbar(lbin, self.fid_bandpowers , yerr=binned_phicl_err , xerr=[lbin - self.lmin, self.lmax - lbin], fmt='o', alpha=0.3)

            plot(ls, phicl, color='k')
            if doN1:
                N1 = self.makeN1(phicl)
                plot(ls, N1 , ls='--', color='r')
            xlim([2, self.phi_lmax])
            # savefig(r'z:\cl_plot.pdf', bbox_inches='tight')

        def chi_squared(self, phicl, clTT=None):

            N1 = self.makeN1(phicl, clTT)
            if clTT is not None:
                phicl = self.makeRenormCPhi(phicl, clTT)
            phicl = phicl + N1 - self.makeN1(lens.fid_phi, None)

            A = np.zeros(self.nbins)
            for b in range(self.nbins):
                A[b] = sum(self.binning_matrix[b, self.lmin[b]:self.lmax[b] + 1] * phicl[self.lmin[b]:self.lmax[b] + 1])
            delta = A - self.bandpowers
            return np.dot(delta, np.dot(self.covinv, delta))


        def testLikes(self):
            print 'chi2 fid=', lens.chi_squared(lens.fid_phi, None)
            ls = self.fid_cl[:2049, 0]
            cl_fid = self.fid_cl[:2049, 1]
            cl2 = cl_fid * (ls / 1000.0) ** 0.03
            phi = self.fid_phi * 0.9
            print 'chi2 =', self.chi_squared(phi, cl2)


if __name__ == "__main__":
    root = 'g60_full_pttptt'

    Duncan = False
    if Duncan:
        f = r'C:\tmp\Planck\lensingApr2014' + os.sep + root + '.dat'
        lens = lensLike(f)
        lens.getRenormDataFromMatrix(r'C:\Work\F90\LensingBiases\like' + os.sep + root + '_renorm_TT_matrix.dat',
                                     r'C:\Work\F90\LensingBiases\like' + os.sep + root + '_renorm_N1_TT_matrix.dat')
        lens.dumpData(r'C:\Work\F90\LensingBiases\like' + os.sep + root)
    else:
        f = r'C:\Work\F90\LensingBiases\like' + os.sep + root + '.dataset'
        lens = lensLike(f)

    lens.testLikes()

    lens.plot()
    show()


if False:
    # check N1s
    ls = np.arange(lens.fid_N1.shape[0])
    norm = (ls * (ls + 1.)) ** 2 / 2 / np.pi
    # plot(ls, lens.fid_N1, color='k')

    x = np.dot(lens.renorm_N1_matrix, lens.fid_cl[lens.renorm_N1_LT, 1])
    plot(lens.renorm_L, x , color='b')
    my = loadtxt('C:\Work\F90\LensingBiases\N1_TT_g60_full_pttptt.dat')
    ls = my[:, 0]

    plot(ls, (ls * (ls + 1.)) ** 2 / 2 / np.pi * my[:, 1], color='r')
    show()
    sys.exit()

if False:
    d = lens.renorm_matrix
    print d.shape
    ls = lens.fid_cl[:2049, 0]
    cl_fid = lens.fid_cl[:2049, 1]
    cl2 = cl_fid * (ls / 1500.0) ** 0.01
    rat = np.zeros(d.shape[0])
    N0 = np.zeros(d.shape[0])
    N0_2 = np.zeros(d.shape[0])

    for i in range(d.shape[0]):
        N0[i] = sum(d[i, 1:] * cl_fid[2:])
        N0_2[i] = sum(d[i, 1:] * cl2[2:])
        rat[i] = (N0[i] / N0_2[i] - 1)
        if (mod(i, 10) == 0): print 'N0, N0_2=', i, sum(d[i, 1:] * cl_fid[2:]), sum(d[i, 1:] * cl2[2:]), rat[i]

    ls = d[:, 0]
    plot(ls, -rat)

