# Clik wrapper from Karim, modified to use standard getdist classes and ClsArray inputs
import clik.smicahlp as smh
import numpy as np
import clik.parobject as php
import clik.hpy as hpy
import clik
import scipy
import os
from CMBlikes import ClsArray
from getdist import ParamNames


class plik_likelihood(object):

    def __init__(self, clikfile, paramnames=None):
        # paramnames will be used to do the mapping to the extra parameters
        self.dffile = clikfile
        self.name = os.path.splitext(os.path.basename(clikfile))[0]
        if isinstance(paramnames, (list, tuple)):
            self.parnames = paramnames
        else:
            if paramnames is None:
                # match versions to the baseline .paramnames file
                for rem in ['_', 'a_', 'b_', 'c_', 'd_']:
                    name = self.name.replace(rem, '_').replace('_bin1', '')
                    paramnames = os.path.join(os.path.dirname(__file__),
                                              '../../data/' + name + '.paramnames')
                    if os.path.exists(paramnames):
                        break
            self.paramnamefile = paramnames
            self.paramnames = ParamNames(paramnames)
            self.parnames = self.paramnames.list()

        self.clik = clik.clik(clikfile)
        self._translate_parname(self.parnames)

        self.fi = hpy.File(self.dffile)

        # some metadata
        self.hascl = self.fi["clik/lkl_0/has_cl"]
        self.lmin = self.fi["clik/lkl_0/lmin"]
        self.lmax = self.fi["clik/lkl_0/lmax"]

        self.mt = self.fi["clik/lkl_0/m_channel_T"] * self.hascl[0]
        self.me = self.fi["clik/lkl_0/m_channel_P"] * self.hascl[1]
        self.mb = self.fi["clik/lkl_0/m_channel_P"] * self.hascl[2]
        self.m = self.mt + self.me + self.mb

        self.nb = self.fi["clik/lkl_0/nbins"] / self.hascl.sum()
        self.rq_shape = (self.nb, self.m, self.m)

        # binning details
        self.blmin = self.fi["clik/lkl_0/bin_lmin"]
        self.blmax = self.fi["clik/lkl_0/bin_lmax"]
        self.b_ws = self.fi["clik/lkl_0/bin_ws"]
        # the binning matrix is also simply obtained this way (but using it is slower than using the binning details, 'cause it's full of zeros)
        self.bns = php.read_bins(self.fi["clik/lkl_0"])

        # compute the binned ells
        self.lm = np.dot(self.bns[:self.nb, :self.lmax - self.lmin + 1], np.arange(self.lmin, self.lmax + 1))

        # get the calibration part (and beam for plik 2015)
        # cal and bal are functions that expect a vector of parameters whose name and ordering are given by cal.varpar and bal.varpar
        # overal calibration is given by cal(pars)*vec(pars)*outer(acmb,acmb)[nm.newaxis,:,:]
        self.cal = smh.calTP_from_smica(self.dffile)
        self.bal = smh.beamTP_from_smica(self.dffile)
        self.acmb = self.fi["clik/lkl_0/A_cmb"]

        # get the binned Cl data array
        self.rqh = self.fi["clik/lkl_0/Rq_hat"]
        self.rqh.shape = self.rq_shape

        # get the additive nuisance components
        self.prms = smh.parametric_from_smica(self.dffile)
        self.prms_name = [p.get_name() for p in self.prms]

        # get the selection vector
        self.oo, self.Jt = smh.ordering_from_smica(self.dffile)

        # get the inverse covariance
        self.siginv = self.fi["clik/lkl_0/criterion_gauss_mat"]
        self.siginv.shape = (len(self.oo), len(self.oo))

        ls = np.arange(self.lmax + 1)
        self.llp1 = ls * (ls + 1) / (2 * np.pi)
        self.llp1[0] = 1
        self.indices = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        self.spectra = ["tt", "ee", "bb", "te", "tb", "eb"]

    def _translate_parname(self, parnames):
        self._translate = {}
        for old, new in zip(self.clik.extra_parameter_names, parnames):
            self._translate[old] = new

    def get_nuisance_vector(self, nuisance_dict, clik_names=None):
        # nuisance dict is a python dict whose keys are defined by the cosmomc paramnames
        if clik_names is None:
            clik_names = self.clik.extra_parameter_names
        return np.array([nuisance_dict[self._translate[par]] for par in clik_names])

    def get_cl_clik_ordering(self, cl_array):
        mcls = []
        lmax = self.lmax
        for i, index in enumerate(self.indices[:4]):
            if self.hascl[i]:
                mcls += [cl_array.get(index)[:lmax + 1] / self.llp1]

        for i in range(4, 6):
            if self.hascl[i]:
                lcls = np.zeros(lmax + 1)
                mcls += [lcls]
        return (mcls)

    def get_clik_vector(self, cl_array, nuisance_dict):
        mcls = self.get_cl_clik_ordering(cl_array)
        extra = self.get_nuisance_vector(nuisance_dict)

        return np.concatenate(mcls + [extra])

    def get_unbinned_nuisance_rq(self, nuisance_dict):
        # compute the nuisance part
        oq = []
        for p in self.prms:
            # prepare the input vector for each component
            pvec = self.get_nuisance_vector(nuisance_dict, p.varpar)
            # compute the rq matrix for this component for those parameters
            oq += [p(pvec)]
            # correct shape for T only component when dealing with T+P case
            if oq[-1].shape[1:] != self.rq_shape[1:]:
                bet = np.zeros((oq[-1].shape[0], self.rq_shape[1], self.rq_shape[1]))
                bet[:, :oq[-1].shape[1], :oq[-1].shape[1]] = oq[-1]
                oq[-1] = bet
        oq = np.array(oq)
        return oq

    def get_nuisance_rq(self, nuisance_dict):
        # get the unbinned one
        oq = self.get_unbinned_nuisance_rq(nuisance_dict)
        oqb = np.zeros((len(oq),) + self.rq_shape)
        # bin it
        nb = self.nb
        m = self.m
        mt = self.mt
        me = self.me
        mb = self.mb
        blmin = self.blmin
        blmax = self.blmax
        b_ws = self.b_ws
        lmin = self.lmin
        lmax = self.lmax
        for b in range(nb):
            if oq.shape[0]:
                oqb[:, b] = np.sum(
                    oq[:, blmin[b]:blmax[b] + 1] * b_ws[np.newaxis, blmin[b]:blmax[b] + 1, np.newaxis, np.newaxis], 1)
        return oqb

    def get_cmb_rq(self, cl_array):
        # get the cls
        mcls = self.get_cl_clik_ordering(cl_array)
        cls = np.zeros((6, self.lmax + 1))
        j = 0
        for i in range(6):
            if self.hascl[i]:
                cls[i] = mcls[j]
                j += 1
        rq = np.zeros((self.nb, self.m, self.m))

        nb = self.nb
        m = self.m
        mt = self.mt
        me = self.me
        mb = self.mb
        blmin = self.blmin
        blmax = self.blmax
        b_ws = self.b_ws
        lmin = self.lmin
        lmax = self.lmax
        # bin it (and order it in)
        for b in range(nb):
            if mt:
                rq[b, :mt, :mt] += np.sum(cls[0, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
                if me:
                    rq[b, :mt, mt:mt + me] += np.sum(
                        cls[3, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
                    rq[b, mt:mt + me, :mt] += np.sum(
                        cls[3, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
                if mb:
                    rq[b, :mt, mt + me:mb + mt + me] += np.sum(
                        cls[4, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
                    rq[b, mt + me:mb + mt + me, :mt] += np.sum(
                        cls[4, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
            if me:
                rq[b, mt:mt + me, mt:mt + me] += np.sum(
                    cls[1, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
                if mb:
                    rq[b, mt:mt + me, mt + me:mb + mt + me] += np.sum(
                        cls[5, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
                    rq[b, mt + me:mb + mt + me, mt:mt + me] += np.sum(
                        cls[5, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
            if mb:
                rq[b, mt + me:mt + me + mb, mt + me:mb + mt + me] += np.sum(
                    cls[2, lmin + blmin[b]:lmin + blmax[b] + 1] * b_ws[blmin[b]:blmax[b] + 1])
        return rq

    def get_calib(self, nuisance_dict):
        # returns a matrix
        cvec = self.get_nuisance_vector(nuisance_dict, self.cal.varpar)
        g = self.cal(cvec)

        cvec = self.get_nuisance_vector(nuisance_dict, self.bal.varpar)
        bg = self.bal(cvec)

        g *= bg
        g *= np.outer(self.acmb, self.acmb)[np.newaxis, :, :]
        return g

    def get_model_vector(self, cl_array, nuisance_dict):
        # get the calib
        g = self.get_calib(nuisance_dict)
        rq = self.get_cmb_rq(cl_array)
        oqb = self.get_nuisance_rq(nuisance_dict)

        rqt = g * (rq + np.sum(oqb, 0))

        return rqt

    def chi_squared(self, cl_array, nuisance_dict, clik=True):
        if clik:
            vec = self.get_clik_vector(cl_array, nuisance_dict)
            return -2 * self.clik(vec)
        else:
            rqt = self.get_model_vector(cl_array, nuisance_dict)
            delta_rq = rqt - self.rqh
            delta = delta_rq.flat[self.oo]

            return np.dot(delta, np.dot(self.siginv, delta))

    def coadd_spectra(self, nuisance_dict):
        good = self.Jt.sum(1) != 0
        Jt = self.Jt[good]
        g = self.get_calib(nuisance_dict)
        Jt = Jt * g.flat[self.oo]
        tm = np.concatenate([self.lm for h in self.hascl if h])

        oqb = self.get_nuisance_rq(nuisance_dict)

        Yo = (g * np.sum(oqb, 0)) - self.rqh
        Yo = Yo.flat[self.oo]

        Jt_siginv = np.zeros((Jt.shape))
        # Jt is full of zero, it's quicker to compute the matrix multiplication this way
        w0, w1 = np.where(Jt != 0)
        for ii in range(len(w0)):
            i = w0[ii]
            j = w1[ii]
            Jt_siginv[i] += Jt[i, j] * self.siginv[j]

        Jt_siginv_Yo = np.dot(Jt_siginv, Yo)

        nl = Jt.shape[0]
        Jt_siginv_J = np.zeros((nl, nl))
        for ii in range(len(w0)):
            j = w0[ii]
            i = w1[ii]
            Jt_siginv_J[j] += Jt[j, i] * Jt_siginv[:, i]
        try:
            rVec = -scipy.linalg.solve(Jt_siginv_J, Jt_siginv_Yo, assume_a='pos')
        except:
            rVec = -np.linalg.solve(Jt_siginv_J, Jt_siginv_Yo)

        tVec = np.zeros(len(tm))
        tVec[good] = rVec
        tm.shape = (-1, len(self.lm))
        tVec.shape = tm.shape
        # beware returns inverse variance of the coadded
        return tm, tVec, Jt_siginv_J

    def coadded_TT(self, dict, want_cov=True):
        ls, coadded, covinv = self.coadd_spectra(dict)
        ls = ls[0, :]
        TT = np.zeros(ls[-1] + 1)
        fac = ls * (ls + 1) / (2 * np.pi)
        TT[ls[0]:] = coadded[0, :] * fac
        if want_cov:
            cov = np.linalg.inv(covinv)[:len(ls), :len(ls)]
            for i in range(cov.shape[0]):
                cov[i, :] *= fac
                cov[:, i] *= fac
            return TT, cov
        else:
            return TT


if __name__ == "__main__":
    import time
    from getdist.types import BestFit

    tag = 'TTTEEE'
    plik = plik_likelihood(
        '/home/aml1005/git/2017/cosmomcplanck/data/clik_14.0/hi_l/plik/plik_rd12_HM_v22b_%s.clik' % tag)
    root = '/scratch/aml1005/Dec17/base/plikHM_%s_lowl_lowE/base_plikHM_%s_lowl_lowE' % (tag, tag)

    fit = BestFit('%s.minimum' % root, want_fixed=True)
    cls = ClsArray('%s.minimum.theory_cl' % root)
    params = fit.getParamDict()
    start = time.time()
    chi2_clik = plik.chi_squared(cls, params)
    print('Likelihood execution time:', time.time() - start)
    for v in fit.chiSquareds:
        if 'plik' in v[1].name:
            test_chi2 = v[1].chisq
            print('Chi-squared calculated: %s, best-fit %s' % (chi2_clik, test_chi2))
            assert (np.isclose(chi2_clik, test_chi2, 0.01))
            break
    if False:  # super slow with unbinned
        start = time.time()
        chi2_python = plik.chi_squared(cls, params, False)
        print('Likelihood execution time:', time.time() - start)
        print(chi2_clik, chi2_python)
    print('Nuisance model', plik.prms_name)
    start = time.time()
    coadded = plik.coadded_TT(params, want_cov=False)
    print('Coadd execution time:', time.time() - start)
