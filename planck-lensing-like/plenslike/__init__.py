import os
import numpy as np
import ctypes as ct

from _plenslike import *

datadir = os.path.dirname(__file__) + "/data/"
pll     = ct.CDLL( os.path.dirname(__file__) + "/_plenslike.so")

# qest
class qest(ct.Structure):
    _fields_ = [ ("ntrm", ct.c_int),
                 ("lmax", ct.c_int),
                 ("s12L", ct.POINTER( ct.POINTER(ct.c_int) )),
                 ("w12L", ct.POINTER( ct.POINTER( ct.POINTER(ct.c_double) ))) ]

pll.fill_qe_resp.argtypes = [ ct.c_int, ct.POINTER(ct.c_double),
                              ct.POINTER(qest), ct.POINTER(qest),
                              ct.POINTER(ct.c_double), ct.c_int,
                              ct.POINTER(ct.c_double), ct.c_int ]

# mono
class plenslike_dat_mono(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmax",          ct.c_int),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("bl_fid",        ct.POINTER(ct.c_double)),
                 ("fl",            ct.POINTER(ct.c_double)),
                 ("vl_inv",        ct.POINTER(ct.c_double)),
                 ("al_inv",        ct.POINTER(ct.c_double)) ]

pll.load_plenslike_dat_mono.argtypes   = [ ct.POINTER(plenslike_dat_mono), ct.c_char_p]
pll.free_plenslike_dat_mono.argtypes   = [ ct.POINTER(plenslike_dat_mono) ]
pll.fill_qe_plm_resp_plm_mono.argtypes = [ ct.c_int, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                                           ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                                           ct.POINTER(ct.c_double), ct.POINTER(ct.c_double) ]
pll.calc_plenslike_mono.restype        = ct.c_double
pll.calc_plenslike_mono_renorm.restype = ct.c_double

class mono():
    def __init__(self, fname):
        print "plenslike:: loading mono likelihood from ", fname

        self.fname = fname
        self.dat = plenslike_dat_mono()
        pll.load_plenslike_dat_mono( ct.byref(self.dat), fname)

    def __del__(self):
        pll.free_plenslike_dat_mono(self.dat)

    def calc_like(self, clpp):
        assert( len(clpp) >= self.dat.lmax )
        return pll.calc_plenslike_mono( ct.byref(self.dat),
                                        clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_like_renorm(self, clpp, cltt, bl):
        assert( len(clpp) >= self.dat.lmax )
        assert( len(cltt) >= self.dat.lmax )
        assert( len(bl)   >= self.dat.lmax )

        return pll.calc_plenslike_mono_renorm( ct.byref(self.dat),
                                               clpp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               cltt.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                               bl.ctypes.data_as(   ct.POINTER(ct.c_double) ) )

    def calc_like_renorm_cltt(self, clpp, cltt):
        bl = np.array( [self.dat.bl_fid[l] for l in xrange(0, self.dat.lmax+1)] )
        return self.calc_like_renorm(clpp, cltt, bl)

    def calc_clpp_renorm_cltt(self, clpp, cltt):
        resp = np.zeros( len(clpp) )
        
        pll.fill_qe_plm_resp_plm_mono( len(resp)-1, resp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                       self.dat.cltt_fid, self.dat.bl_fid, self.dat.fl,
                                       cltt.ctypes.data_as( ct.POINTER(ct.c_double) ), self.dat.bl_fid )

        return resp**2 / np.array( [ self.dat.al_inv[l] for l in xrange(0, len(clpp)) ] )**2 * clpp

    def calc_bins_clpp(self, clpp):
        bins = np.zeros( self.dat.nbins )

        pll.fill_plenslike_mono_bins( ct.byref(self.dat),
                                      bins.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                      clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

        return bins

# quad
class plenslike_dat_quad(ct.Structure):
    _fields_ = [ ("nbins",         ct.c_int),
                 ("lmax",          ct.c_int),
                 ("lmaxt",         ct.c_int),
                 ("lmax1",         ct.c_int),
                 ("lmax2",         ct.c_int),
                 ("lmax3",         ct.c_int),
                 ("lmax4",         ct.c_int),
                 ("s4hat",         ct.c_double),
                 ("s4std",         ct.c_double),
                 ("bin_lmins",     ct.POINTER(ct.c_int)),
                 ("bin_lmaxs",     ct.POINTER(ct.c_int)),
                 ("bin_vals",      ct.POINTER(ct.c_double)),
                 ("mat_sigma",     ct.POINTER(ct.c_double)),
                 ("mat_sigma_inv", ct.POINTER(ct.c_double)),
                 ("clpp_fid",      ct.POINTER(ct.c_double)),
                 ("vl_inv",        ct.POINTER(ct.c_double)),
                 ("rl_inv",        ct.POINTER(ct.c_double)),
                 ("sl_fid",        ct.POINTER(ct.c_double)),
                 ("cltt_fid",      ct.POINTER(ct.c_double)),
                 ("bl1n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl2n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl3n1_fid",     ct.POINTER(ct.c_double)),
                 ("bl4n1_fid",     ct.POINTER(ct.c_double)),
                 ("qe12",          ct.POINTER(qest)),
                 ("qe34",          ct.POINTER(qest)) ]

pll.load_plenslike_dat_quad.argtypes   = [ ct.POINTER(plenslike_dat_quad), ct.c_char_p]
pll.free_plenslike_dat_quad.argtypes   = [ ct.POINTER(plenslike_dat_quad) ]

pll.fill_quad_resp_pp_cltt.argtypes    = [ ct.c_int, ct.POINTER(ct.c_double), ct.POINTER(plenslike_dat_quad), ct.POINTER(ct.c_double) ]

pll.calc_plenslike_quad.restype        = ct.c_double
pll.calc_plenslike_quad_renorm_cltt.restype = ct.c_double

class quad():
    def __init__(self, fname):
        print "plenslike:: loading quad likelihood from ", fname

        self.fname = fname
        self.dat = plenslike_dat_quad()
        pll.load_plenslike_dat_quad( ct.byref(self.dat), fname)

    def __del__(self):
        pll.free_plenslike_dat_quad(self.dat)

    def calc_qc_resp_pp_cltt(self, lmax, cltt):
        assert( len(cltt) >= self.dat.lmaxt )
        
        ret = np.zeros(lmax+1)
        pll.fill_quad_resp_pp_cltt( lmax, ret.ctypes.data_as(ct.POINTER(ct.c_double)),
                                    ct.byref(self.dat), cltt.ctypes.data_as(ct.POINTER(ct.c_double)) )
        return ret

    def calc_like(self, clpp):
        assert( len(clpp) >= self.dat.lmax )
        return pll.calc_plenslike_quad( ct.byref(self.dat),
                                        clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_like_renorm_cltt(self, clpp, cltt):
        assert( len(clpp) >= self.dat.lmax )
        assert( len(cltt) >= self.dat.lmaxt )

        return pll.calc_plenslike_quad_renorm_cltt( ct.byref(self.dat),
                                                    clpp.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                                    cltt.ctypes.data_as( ct.POINTER(ct.c_double) ) )

    def calc_clpp_renorm_cltt(self, clpp, cltt):
        resp = self.calc_qc_resp_pp_cltt( len(clpp)-1, cltt )
        return resp * np.array( [ self.dat.rl_inv[l] for l in xrange(0, len(clpp)) ] ) * clpp

    def calc_bins_clpp(self, clpp):
        bins = np.zeros( self.dat.nbins )

        pll.fill_plenslike_quad_bins( ct.byref(self.dat),
                                      bins.ctypes.data_as( ct.POINTER(ct.c_double) ),
                                      clpp.ctypes.data_as( ct.POINTER(ct.c_double) ) )

        return bins
