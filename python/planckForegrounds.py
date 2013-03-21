import paramNames
import numpy as np


class foregroundModel():
    def __init__(self, freq_params, lmax):

        self.lmax = lmax
        self.ls = range(lmax + 1)
        self.A_ps_100 = freq_params[0]
        self.A_ps_143 = freq_params[1]
        self.A_ps_217 = freq_params[2]
        self.A_cib_143 = freq_params[3]
        self.A_cib_217 = freq_params[4]
        self.A_sz = freq_params[5]
        self.r_ps = freq_params[6]
        self.r_cib = freq_params[7]
        self.ncib = freq_params[8]
        self.cal0 = freq_params[9]
        self.cal1 = freq_params[10]
        self.cal2 = freq_params[11]
        self.xi = freq_params[12]
        self.A_ksz = freq_params[13]

    def CIB_cl(self):
        cl_cib = (self.ls / 3000) ** (self.ncib)
        return cl_cib


