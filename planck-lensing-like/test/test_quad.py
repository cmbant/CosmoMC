#!/usr/bin/env python
# simple script to test loading and likelihood calculation for 'mono' format files.

import glob
import numpy as np
import pylab as pl
import scipy.stats

import plenslike

for like in [ plenslike.quad(tfname) for tfname in glob.glob(plenslike.datadir + "dx9nom/*quad*.dat") ]:
    dat = like.dat

    cltt   = np.array( [like.dat.cltt_fid[l] for l in xrange(0, like.dat.lmaxt+1) ] )
    clpp   = np.array( [like.dat.clpp_fid[l] for l in xrange(0, like.dat.lmax+1) ] )
    rln1   = np.array( [like.dat.rl_inv[l] for l in xrange(0, like.dat.lmax+1) ] )

    print "tfname = ", like.fname.replace(plenslike.datadir, '')
    print "   nbins = ", dat.nbins
    print "   lmax  = ", dat.lmax, " lmaxt = ", dat.lmaxt
    print "      s4 = ", dat.s4hat, " +- ", dat.s4std
    print "   lmins = ", [dat.bin_lmins[i] for i in xrange(0,dat.nbins)]
    print "   lmaxs = ", [dat.bin_lmaxs[i] for i in xrange(0,dat.nbins)]
    print "   vals  = ", [dat.bin_vals[i]  for i in xrange(0,dat.nbins)]
    print "   ntrm  = ", dat.qe12[0].ntrm, " (qe12) ", dat.qe34[0].ntrm, " (qe34)"
    print "   PTE   = ", 100.*(1.-scipy.stats.chi2.cdf( like.calc_like(clpp)*-2., like.dat.nbins))

    rln1_recalc = 1./like.calc_qc_resp_pp_cltt( dat.lmax, cltt )

    ls = np.arange(0, dat.lmax+1)
    pl.loglog((ls*(ls+1.))**2 / (2.*np.pi) * np.sqrt(rln1[0:dat.lmax+1]))
    pl.loglog((ls*(ls+1.))**2 / (2.*np.pi) * np.sqrt(rln1_recalc[0:dat.lmax+1]), ls='--', lw=2)

pl.ion()
pl.show()
