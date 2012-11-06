import os
import ctypes as ct

from _plenslike import *

datadir = os.path.dirname(__file__) + "/data/"

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

pll = ct.CDLL( os.path.dirname(__file__) + "/_plenslike.so")
pll.load_plenslike_dat_mono.argtypes   = [ ct.POINTER(plenslike_dat_mono), ct.c_char_p]
pll.free_plenslike_dat_mono.argtypes   = [ ct.POINTER(plenslike_dat_mono) ]
pll.calc_plenslike_mono.restype        = ct.c_double
pll.calc_plenslike_mono_renorm.restype = ct.c_double

class mono():
    def __init__(self, tfname):
        print "plenslike:: loading mono likelihood from ", tfname

        self.dat = plenslike_dat_mono()
        pll.load_plenslike_dat_mono( ct.byref(self.dat), tfname)

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
