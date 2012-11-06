#!/usr/bin/env python
# simple script to test loading and likelihood calculation for 'mono' format files.

import numpy as np
import plenslike

for mono in [ plenslike.mono(plenslike.datadir + "v50/plenslike_v50_mono_143.dat"),
              plenslike.mono(plenslike.datadir + "v50/plenslike_v50_mono_tot_143.dat") ]:

    dat = mono.dat

    print "nbins = ", dat.nbins
    print "lmax  = ", dat.lmax
    print "lmins = ", [dat.bin_lmins[i] for i in xrange(0,dat.nbins)]
    print "lmaxs = ", [dat.bin_lmaxs[i] for i in xrange(0,dat.nbins)]
    print "vals  = ", [dat.bin_vals[i]  for i in xrange(0,dat.nbins)]

    print "sigma = "
    for i in xrange(0, dat.nbins):
        print [dat.mat_sigma[i*dat.nbins + j] for j in xrange(0,dat.nbins)]

    print "sigma_inv = "
    for i in xrange(0, dat.nbins):
        print [dat.mat_sigma_inv[i*dat.nbins + j] for j in xrange(0,dat.nbins)]

    lmax = dat.lmax
    clpp_fid = np.array( [dat.clpp_fid[l] for l in xrange(0,lmax+1)] )
    cltt_fid = np.array( [dat.cltt_fid[l] for l in xrange(0,lmax+1)] )
    bl_fid   = np.array( [dat.bl_fid[l] for l in xrange(0,lmax+1)] )
    fl       = np.array( [dat.fl[l] for l in xrange(0,lmax+1)] )
    vl_inv   = np.array( [dat.vl_inv[l] for l in xrange(0,lmax+1)] )
    al_inv   = np.array( [dat.al_inv[l] for l in xrange(0,lmax+1)] )

    print "like = ", mono.calc_like(clpp_fid)
    print "like renorm = ", mono.calc_like_renorm(clpp_fid, cltt_fid, bl_fid)
