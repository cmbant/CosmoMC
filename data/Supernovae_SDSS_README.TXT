                 July 2009

Enclosed are the light curve fit results corresponding
to the (Kessler et al 2009) ApJS paper:

  "First-year Sloan Digital Sky Survey-II (SDSS-II) Supernova 
   Results: Hubble Diagram and Cosmological Parameters"


The MLCS.FITRES file contains the fit results ordered by
sample (LOWZ, SDSS, ESSENCE, SNLS, HST). The header is 
self-documented so that it can be easiliy parsed.  
The sample-combinations in the paper can be obtained
by simply combining the appropriate SNe.

The SALT2.FITRES file contains the analogous SALT2 results as a 
function of sample. However, in contrast to MLCS, the distance 
modulus is given for each sample combination, MU[a-f] and MUERR[a-f].
The distance & error are set to -9 for SNe that are not included
in the particular sample-combination; for example, MUa = -9 for
all but the SDSS SN. As a reminder, the sample-combinations are

 a : SDSS only
 b : SDSS+ESSENCE+SNLS 
 c : LOWZ+SDSS 
 d : LOWZ+SDSS+ESSENCE+SNLS 
 e : LOWZ+SDSS+ESSENCE+SNLS+HST 
 f : LOWZ+ESSENCE+SNLS 


The column header-words for MLCS.FITRES are defined below:

CID      : candidate id, or SN id (character or integer)
Z        : redshift
ZERR     : error on redshift
MU       : distance modulus
MUERR    : error on distance modulus
AV       : host galaxy extinction in V-band
AVERR    : error on above
DELTA    : MLCS shape-luminosity parameter
DELTAERR : error on above
PKMJD    : MJD at maximum brightness of rest-frame B-band
PKMJDERR : error on above
RV       : CCM89 parameter for extinction law
RVERR    : error on above
CHI2     : mlcs fit-chi2
NDOF     : Number of degrees of freedom in fit
IDTEL    : telescope ID (1=SDSS, 3=ESSE, 4=SNLS, 50=Nearby, 100=HST)


and for SALT2.FITRES we have
c     : color term
cerr  : error on above
x1    : shape-luminosoty parameter
x1err : error on above

=========== END ===========
