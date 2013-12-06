#Make BBN interpolation table from AlterBBN
#Used AlterBBN 1.4, compiled after changing exit code to 0, and increasing sf of the text output in alter_etannutau

import numpy as np
import subprocess
import sys

hbar = 1.05457e-34
c = 299792458.
kB = 1.380650e-23
MeV = 1.6021764e-13
eV = MeV / 1e6
G = 6.6738e-11
TCMB = 2.7255
mP = np.sqrt(hbar * c / G)
Mpc = 3.0856776e22

m_proton = 1.672621637e-27
m_H = 1.673575e-27
not4 = 3.9715
m_He = m_H * not4

zeta3 = 1.202056903
tau_n = 880.1

n_photon = (kB * TCMB / hbar / c) ** 3 * zeta3 * 2 / np.pi ** 2

omegafac = (1e5 / Mpc) ** 2 / (8 * np.pi * G) * 3


def runBBN(eta,DeltaN=0,tau=tau_n,prog='./alter_etannutau.x'):
    return subprocess.check_output([prog, str(eta), str(3+DeltaN), str(tau)], shell=False,stderr=subprocess.STDOUT)

    
def runBBNplot():
    #standard BBN plot
    etas = 10**np.linspace(-12,-8, num=50, endpoint=True)
    lines=[]
    for eta in etas:
        outtxt = runBBN(eta,prog='./main.x')
        (YBBN, D, He3, Li7, Li6, Be7) = [float(s.strip()) for s in outtxt.split('\n')[0].split()]
        line = ('%12.5e'*7)%(eta, YBBN, D, He3, Li7, Li6, Be7) 
        print line
        lines += [line]
    f=open('bbn_abundance_eta.dat', 'w')
    f.write("\n".join(lines))
    f.close()
    sys.exit()

        
ombh2s = np.hstack((np.linspace(0.005, 0.02, num=int((0.02 - 0.005) / 0.001), endpoint=False),
         np.linspace(0.020,0.024, num=16 , endpoint=False),
         np.linspace(0.024, 0.04, num=int((0.04 - 0.024) / 0.001) + 1, endpoint=True)))

DeltaNs = [-3,-2,-1,-0.5,0,0.5,1,2,3,4,5,6,7]

lines=[]
for DeltaN in DeltaNs:
# this is fiducial mass fraction. For low ombh2 has no effect, for subsequent steps use previous value
    Yp_fid = 0.2;
    for ombh2 in ombh2s:
        rho_b = ombh2 * omegafac;
        n_baryon = (4 * Yp_fid / m_He + (1 - Yp_fid) / m_H) * rho_b
        eta = n_baryon / n_photon


        outtxt = runBBN(eta,DeltaN)
        
        (YBBN, D, He3, Li7, Li6, Be7) = [float(s.strip()) for s in outtxt.split('\n')[7].split()[1:]]
        (sigYBBN, sigdD, sigHe3, sigLi7, sigLi6, sigBe7) = [float(s.strip()) for s in outtxt.split('\n')[8].split()[2:]]

        actual_ombh2 = eta * n_photon * (YBBN * m_He / 4 + (1 - YBBN) * m_H)/omegafac
        Yp= (eta * n_photon/ actual_ombh2/omegafac- 1/m_H)/(4/m_He-1/m_H)
        Yp_fid = Yp
#        print ombh2, DeltaN, eta, YBBN, actual_ombh2, Yp
        line = (('%12.5f ')*6 + ('%12.3e %12.2e')*3)%(actual_ombh2, eta*1e10, DeltaN, Yp, YBBN, sigYBBN, D, sigdD, He3, sigHe3,Li7,sigLi7) 
        lines += [line]
        print line
    lines+=['']

f=open('BBN_full_alterBBN_'+str(tau_n)+'.dat', 'w')
f.write('''#BBN prediction of the primordial Helium abundance $Y_p$ as 
#function of the baryon density $\omega_b h^2$ and number of 
#extra radiation degrees of freedom $\Delta N$.
#Calculated with AlterBBN v1.4 [http://superiso.in2p3.fr/relic/alterbbn/] for a 
#neutron lifetime of %.1f s and CMB temperature %.4f K
#Yp^BBN is the BBN-standard nucleon number fraction, Yp is the mass fraction for CMB codes

#      ombh2        eta10       DeltaN           Yp       Yp^BBN     sigma_Yp          D/H     err D/H       He3/H    err He3/H         Li7      sig Li7
 
'''%(tau_n,TCMB))
f.write("\n".join(lines))
f.close()

