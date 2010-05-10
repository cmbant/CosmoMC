#Settings here are for Intel ifort on CITA cluster; uncomment others below
#You may need to edit the library paths for MKL for Intel
#Beware of using optmizations which lose accuracy - may give errors when running
#With Visual Fortan open cosmomc.dsw workspace file instead.

#Edit for CFITSIO, HEALPIX directories
cfitsio = /home/pascal/include/cfitsio
healpix = /home/pascal/include/Healpix_2.00
WMAP3 = /home/antlewis/WMAP3

IFLAG = -I

#Intel MPI (assuming mpif77 set to point to ifc)
#change -lmkl_ia32 to -lmkl_p3 for MKL versions earlier than 6 (6 needed for ifc 8+) 
F90C     = mpif77
FFLAGS = -O2 -Vaxlib -ip -W0 -WB -openmp -fpp -DMPI
LAPACKL = -L/opt/intel/mkl61/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread

# Use with Intel 9 on the Altix (eg. COSMOS)
# (do "module load icomp90" before compiling, or add to .bashrc)
#F90C = ifort
#FFLAGS = -O3 -w -fpp2 -DMPI 
#LAPACKL = -lscs -lguide -lmpi -ldl

#Intel fortran 8, check you have the latest update from the Intel web pages
#See Makefile_intel for ifc 7.1 or lower (some versions have problems)
#F90C     = ifort
#FFLAGS = -O2 -Vaxlib -ip -W0 -WB -openmp -fpp
#LAPACKL = -L/opt/intel/mkl61/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread

#Intel fortran 8 on COSMOS
#F90C     = ifort
#FFLAGS = -O2 -Vaxlib -W0 -WB -openmp -fpp2 -lmpi -lscs  -ldl -DMPI
#LAPACKL = -L/opt/intel/mkl80/lib/64 -lmkl_lapack -lmkl_ipf -limf

#SGI, -mp toggles multi-processor. Use -O2 if -Ofast gives problems.
#Not various versions of the compiler are buggy giving erroneous seg faults with -mp.
#Version 7.3 is OK, currently version 7.4 is bugged, as are some earlier versions.
#F90C     = f90
#LAPACKL = -lcomplib.sgimath
#FFLAGS  = -Ofast -mp

#Digital/Compaq fortran, -omp toggles multi-processor
#F90C    = f90
#FFLAGS  = -omp -O4 -arch host -math_library fast -tune host -fpe1
#LAPACKL  = -lcxml

#Intel Itanium efc, -openmp toggles multi-processor:
#F90C     = efc
#FFLAGS = -O2 -Vaxlib -W0 -WB -openmp -fpp
#LAPACKL = -L/opt/intel/mkl/lib/64 -lmkl_lapack -lmkl_ipf -lguide -lpthread -static

#Intel Itanium efc MPI:
# mpirun -np 4 ./cosmomc params.ini
#F90C     = efc
#FFLAGS = -O2 -Vaxlib -W0 -WB -openmp -fpp -lmpi -lscs  -ldl -DMPI
#LAPACKL = -L/opt/intel/mkl/lib/64 -lmkl_lapack -lmkl_ipf -limf

#G95; make sure LAPACK and MPI libs also compiled with g95
#F90C     = mpif90
#FFLAGS = -O2 -DMPI
#LAPACKL = /LAPACK/lapack_LINUX.a /LAPACK/blas_LINUX.a

#Absoft ProFortran, single processor, set -cpu:[type] for your local system 
#F90C     = f95
#FFLAGS = -O2 -s -cpu:athlon -lU77 -w -YEXT_NAMES="LCS" -YEXT_SFX="_" 
#LAPACKL =  -llapack -lblas -lg2c 
#IFLAG = -p

#NAGF95, single processor:
#F90C     = f95
#FFLAGS = -DNAGF95 -O3
#LAPACKL = -llapack -lblas -lg2c

#Sun, single processor:
#F90C     = f90
#FFLAGS = -fast -ftrap=%none
#LAPACKL = -lsunperf -lfsu
#LAPACKL = -lsunperf -lfsu -lsocket -lm
#IFLAG = -M

#Sun MPI
#F90C = mpf90
#FFLAGS =  -O4 -openmp -ftrap=%none -dalign -DMPI
#LAPACKL = -lsunperf -lfsu -lmpi_mt
#IFLAG = -M

#Sun parallel enterprise:
#F90C     = f95
#FFLAGS =  -O4 -xarch=native64 -openmp -ftrap=%none
#LAPACKL = -lsunperf -lfsu
#IFLAG = -M


#IBM XL Fortran, multi-processor (run "module load lapack" then run "gmake")
# See also http://cosmocoffee.info/viewtopic.php?t=326
#F90C     = xlf90_r $(LAPACK)
#FFLAGS  = -WF,-DIBMXL -qsmp=omp -qsuffix=f=f90:cpp=F90 -O3 -qstrict -qarch=pwr3 -qtune=pwr3
#INCLUDE = -lessl
#LAPACKL =

PROPOSE = propose.o
CLSFILE = CMB_Cls_simple.o

#Can use params_H if you prefer more generic parameters
PARAMETERIZATION = params_CMB.o

WMAP3_inc = $(IFLAG)$(cfitsio)/include $(IFLAG)$(healpix)/include $(IFLAG)$(WMAP3)
WMAP3_lib = -L$(cfitsio)/lib -L$(healpix)/lib -L$(WMAP3) \
	-lhealpix -lcfitsio

CLSLIB  = $(WMAP3_lib) -L../camb -lcamb $(LAPACKL) 

F90FLAGS = $(FFLAGS) $(WMAP3_inc) $(IFLAG)../camb $(INCLUDE)

DISTFILES = utils.o settings.o GetDist.o

WMAPOBJS = $(WMAP3)/read_archive_map.o \
	$(WMAP3)/read_fits.o \
	$(WMAP3)/WMAP_3yr_options.o \
	$(WMAP3)/WMAP_3yr_util.o \
	$(WMAP3)/WMAP_3yr_tt_pixlike.o \
	$(WMAP3)/WMAP_3yr_tt_beam_and_ptsrc_corr.o \
	$(WMAP3)/WMAP_3yr_teeebb_pixlike.o \
	$(WMAP3)/WMAP_3yr_likelihood.o

OBJFILES= $(WMAPOBJS) utils.o settings.o cmbtypes.o cmbdata.o WeakLen.o mpk.o \
	supernovae.o SDSSLy-a-v3.o\
	$(CLSFILE) paramdef.o $(PROPOSE) $(PARAMETERIZATION) calclike.o \
	conjgrad_wrapper.o EstCovmat.o postprocess.o MCMC.o driver.o

default: cosmomc

all : cosmomc getdist

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


cosmomc: $(OBJFILES) 	
	$(F90C) -o ../cosmomc $(OBJFILES) $(CLSLIB) $(F90FLAGS)


clean:
	rm -f *.o *.mod *.d *.pc *.obj ../core

getdist: $(DISTFILES)
	$(F90C) -o ../getdist $(DISTFILES) $(CLSLIB) $(F90FLAGS) 
