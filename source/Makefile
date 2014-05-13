
OUTPUT_DIR ?= Release

BUILD ?= MPI
#set BUILD to MPI to force MPI, should be set in ../Makefile

ifortErr = $(shell which ifort >/dev/null; echo $$?)

#these settings for ifort 13 and higher
ifeq "$(ifortErr)" "0"
#ifort; Can remove -xHost if your cluster is not uniform, or specify specific processor optimizations -x...
F90C     = ifort
MPIF90C = mpif90
FFLAGS = -mkl -openmp -O3 -xHost -no-prec-div -fpp
DEBUGFLAGS = -mkl -openmp -g -check all -check noarg_temp_created -traceback -no-prec-div -fpp 
#add -fpe0 to check for floating point errors (think lowLike also throws these harmlessly)
MODOUT = -module $(OUTPUT_DIR)
LAPACKL =
else
#gfortran does not yet work yet, open GCC 4.9 compiler bugs are
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60414
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60334
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60289
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=60255
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=51976
#and preferably
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=44672
F90C = gfortran
MPIF90C = mpif90
FFLAGS = -cpp -O3 -ffast-math
DEBUGFLAGS = -cpp -Og -g -fbounds-check -fbacktrace
MODOUT =  -J$(OUTPUT_DIR)
#LAPACKL = -Wl,-framework -Wl,accelerate
#commented above is (I think) for Mac; this is standard linux (sudo apt-get install liblapack-dev)
LAPACKL = -lblas -llapack
endif

#use "make RECOMBINATION=cosmorec" to build with CosmoRec rather than RECFAST default
RECOMBINATION ?=recfast

EXTDATA=

#use PRECISION=SINGLE to use single precision (probably doesn't work)
PRECISION ?=

#set WMAP empty not to compile with WMAP, e.g. WMAP = /scratch/../WMAP9/likelihood_v5
WMAP ?=

#set PICO=YES to build with pico support, or set path to PICO manually
#e.g. PICO=/path/to/pypico
PICO ?=
ifeq ($(PICO),YES)
override PICO =$(shell python -c "import pypico; print pypico.get_folder()")
endif


#Needed for WMAP; set from ~/.bashrc definition or edit here
cfitsio ?=

#empty not to compile with CLIK, set from ~/.bashrc
PLANCKLIKE ?=

#These is not used in public cosmomc
NONCLIKLIKE ?=
highL ?=
#highL = ../highL

IFLAG = -I
INCLUDE =

#clik path and library
ifeq ($(PLANCKLIKE),cliklike)
CLIKL = -L$(CLIKPATH)/lib -lclik_f90
INCLUDE = $(IFLAG)$(CLIKPATH)/include
endif


#COSMOS: use "module load cosmolib latest"
#use "runCosmomc" (globally installed) to run, defining required memory usage
ifeq ($(COSMOHOST),cosmos)
F90C = ifort
FFLAGS = -openmp -fast -w -fpp2
LINKMPI = -lmpi
LAPACKL = -mkl=sequential -lmkl_lapack
cfitsio = $(CFITSIO)
WMAP = $(COSMOLIB)/WMAP9
GSLPATH = $(GSL_ROOT)
endif

ifeq ($(COSMOHOST),darwin)
WMAP = /scratch/aml1005/wmap_likelihood_v5
cfitsio = /usr/local/Cluster-Apps/cfitsio/intel/3.300
endif

ifeq ($(OUTPUT_DIR),Debug)
FFLAGS=$(DEBUGFLAGS)
endif


ifeq ($(BUILD),MPI)
override OUTPUT_DIR :=$(OUTPUT_DIR)MPI
FFLAGS +=  -DMPI
LAPACKL +=  $(LINKMPI)
F90C = $(MPIF90C)
endif

#would like to embed this somehow..
#GIT_HASHTAG = "git log --pretty=format:'%H:%cd' -n 1"

CALCFILES = $(OUTPUT_DIR)/Calculator_CAMB.o

#Can use params_H if you prefer more generic parameters
PARAMETERIZATION = $(OUTPUT_DIR)/CosmologyParameterizations.o

GETDISTFLAGS = $(FFLAGS)
F90FLAGS = $(FFLAGS) $(IFLAG)../camb $(INCLUDE)
LINKFLAGS = -L../camb -lcamb_$(RECOMBINATION) $(LAPACKL) $(F90CRLINK) $(CLIKL)

ifneq ($(PICO),)
CALCFILES += $(OUTPUT_DIR)/fpico_interface.o $(OUTPUT_DIR)/Calculator_PICO.o
LINKFLAGS +=  $(shell python -c "import pypico; print pypico.get_link()")
F90FLAGS += -DPICO $(shell python -c "import pypico; print pypico.get_include()")
endif


UTILS = $(OUTPUT_DIR)/MiscUtils.o $(OUTPUT_DIR)/ArrayUtils.o $(OUTPUT_DIR)/StringUtils.o $(OUTPUT_DIR)/MpiUtils.o $(OUTPUT_DIR)/FileUtils.o \
    $(OUTPUT_DIR)/RandUtils.o $(OUTPUT_DIR)/ObjectLists.o $(OUTPUT_DIR)/Interpolation.o

DISTFILES = $(UTILS) $(OUTPUT_DIR)/IniObjects.o \
    $(OUTPUT_DIR)/ObjectParamNames.o $(OUTPUT_DIR)/Matrix_utils_new.o \
    $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/samples.o $(OUTPUT_DIR)/IO.o $(OUTPUT_DIR)/GetDist.o


ifneq ($(PRECISION),)
FFLAGS += -D$(PRECISION) -DMATRIX_$(PRECISION)
endif


ifneq ($(PLANCKLIKE),)
FFLAGS += -DCLIK
endif

ifneq ($(highL),)
FFLAGS += -DhighL
LINKFLAGS += -L$(highL) -lhigh
INCLUDE += $(IFLAG)$(highL)
endif

PLANCKLIKEFILES=


ifneq ($(PLANCKLIKE),)
PLANCKLIKEFILES += $(OUTPUT_DIR)/$(PLANCKLIKE).o
endif

LIKEFILES =  $(OUTPUT_DIR)/DataLikelihoods.o $(OUTPUT_DIR)/calclike.o

OBJFILES = $(UTILS) $(OUTPUT_DIR)/IniObjects.o $(OUTPUT_DIR)/ObjectParamNames.o \
	$(OUTPUT_DIR)/Matrix_utils_new.o  $(OUTPUT_DIR)/settings.o \
	$(OUTPUT_DIR)/samples.o $(OUTPUT_DIR)/IO.o $(OUTPUT_DIR)/likelihood.o $(OUTPUT_DIR)/GeneralTypes.o $(OUTPUT_DIR)/BaseParameters.o \
	$(OUTPUT_DIR)/propose.o $(OUTPUT_DIR)/ParamSet.o $(OUTPUT_DIR)/PowellConstrainedMinimize.o $(OUTPUT_DIR)/ImportanceSampling.o \
	$(OUTPUT_DIR)/MCMC.o $(OUTPUT_DIR)/minimize.o $(OUTPUT_DIR)/SampleCollector.o $(OUTPUT_DIR)/GeneralSetup.o


SUPERNOVAE = $(OUTPUT_DIR)/supernovae_Union2.o $(OUTPUT_DIR)/supernovae_SNLS.o $(OUTPUT_DIR)/supernovae_JLA.o

DATAMODULES = $(PLANCKLIKEFILES) $(OUTPUT_DIR)/mpk.o $(OUTPUT_DIR)/wigglez.o \
	$(OUTPUT_DIR)/bao.o $(SUPERNOVAE) $(OUTPUT_DIR)/supernovae.o $(OUTPUT_DIR)/HST.o  $(OUTPUT_DIR)/CMB.o $(OUTPUT_DIR)/CMBlikes.o


COSMOFILES = $(OUTPUT_DIR)/CosmologyTypes.o $(OUTPUT_DIR)/CosmoTheory.o $(OUTPUT_DIR)/bbn.o $(OUTPUT_DIR)/Calculator_Cosmology.o \
	$(OUTPUT_DIR)/CalcLike_Cosmology.o $(OUTPUT_DIR)/Likelihood_Cosmology.o $(CALCFILES) $(DATAMODULES) $(PARAMETERIZATION) $(LIKEFILES) \
	$(OUTPUT_DIR)/CosmologyConfig.o

OBJFILES += $(COSMOFILES) $(OUTPUT_DIR)/driver.o

F90CRLINK =

ifeq ($(RECOMBINATION),cosmorec)
## This is flag is passed to the Fortran compiler allowing it to link C++ (uncomment the right one).
# GCC (gfortran/g++)
COSMOREC_PATH ?= ../CosmoRec/
F90CRLINK = -L$(COSMOREC_PATH) -lCosmoRec -L$(GSLPATH)/lib -lgsl -lgslcblas -lstdc++
# Intel Compilers (ifort/icpc)
#F90CRLINK = -cxxlib -L$(COSMOREC_PATH) -lCosmoRec -L$(GSLPATH)/lib -lgsl -lgslcblas
FFLAGS +=  -DCOSMOREC
endif

ifeq ($(RECOMBINATION),hyrec)
HYREC_PATH ?= ../HyRec/
F90CRLINK += -L$(HYREC_PATH) -lhyrec
endif

default: cosmomc

$(OUTPUT_DIR)/RandUtils.o: $(OUTPUT_DIR)/MpiUtils.o
$(OUTPUT_DIR)/FileUtils.o: $(OUTPUT_DIR)/MpiUtils.o $(OUTPUT_DIR)/MiscUtils.o $(OUTPUT_DIR)/StringUtils.o
$(OUTPUT_DIR)/ObjectLists.o: $(OUTPUT_DIR)/FileUtils.o
$(OUTPUT_DIR)/IniObjects.o: $(OUTPUT_DIR)/FileUtils.o
$(OUTPUT_DIR)/ObjectParamNames.o: $(UTILS)
$(OUTPUT_DIR)/settings.o:  $(OUTPUT_DIR)/IniObjects.o $(OUTPUT_DIR)/ObjectParamNames.o
$(OUTPUT_DIR)/propose.o: $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/Matrix_utils_new.o
$(OUTPUT_DIR)/IO.o: $(OUTPUT_DIR)/settings.o
$(OUTPUT_DIR)/samples.o: $(OUTPUT_DIR)/ObjectLists.o $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/Matrix_utils_new.o
$(OUTPUT_DIR)/GetDist.o: $(OUTPUT_DIR)/IO.o $(OUTPUT_DIR)/samples.o

$(OUTPUT_DIR)/likelihood.o: $(OUTPUT_DIR)/GeneralTypes.o
$(OUTPUT_DIR)/GeneralTypes.o: $(OUTPUT_DIR)/ObjectLists.o $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/IO.o
$(OUTPUT_DIR)/BaseParameters: $(OUTPUT_DIR)/GeneralTypes.o $(OUTPUT_DIR)/IO.o
$(OUTPUT_DIR)/ParamSet.o: $(OUTPUT_DIR)/BaseParameters.o $(OUTPUT_DIR)/samples.o
$(OUTPUT_DIR)/DataLikelihoods.o: $(OUTPUT_DIR)/likelihood.o $(OUTPUT_DIR)/ParamSet.o $(DATAMODULES)
$(OUTPUT_DIR)/calclike.o: $(OUTPUT_DIR)/DataLikelihoods.o
$(OUTPUT_DIR)/ImportanceSampling.o: $(OUTPUT_DIR)/calclike.o $(OUTPUT_DIR)/IO.o
$(OUTPUT_DIR)/MCMC.o: $(OUTPUT_DIR)/calclike.o $(OUTPUT_DIR)/RandUtils.o  $(OUTPUT_DIR)/propose.o $(OUTPUT_DIR)/ParamSet.o
$(OUTPUT_DIR)/minimize.o: $(OUTPUT_DIR)/PowellConstrainedMinimize.o $(OUTPUT_DIR)/MCMC.o
$(OUTPUT_DIR)/SampleCollector.o: $(OUTPUT_DIR)/MCMC.o
$(OUTPUT_DIR)/GeneralSetup.o: $(OUTPUT_DIR)/SampleCollector.o $(OUTPUT_DIR)/ImportanceSampling.o  $(OUTPUT_DIR)/minimize.o

$(OUTPUT_DIR)/driver.o: $(OUTPUT_DIR)/GeneralSetup.o $(OUTPUT_DIR)/CosmologyConfig.o

$(PARAMETERIZATION): $(OUTPUT_DIR)/ParamSet.o $(OUTPUT_DIR)/bbn.o
$(OUTPUT_DIR)/Calculator_Cosmology.o: $(OUTPUT_DIR)/CosmoTheory.o
$(OUTPUT_DIR)/CosmoTheory.o: $(OUTPUT_DIR)/CosmologyTypes.o $(OUTPUT_DIR)/likelihood.o
$(OUTPUT_DIR)/supernovae.o:  $(OUTPUT_DIR)/Likelihood_Cosmology.o $(SUPERNOVAE)
$(OUTPUT_DIR)/Calculator_CAMB.o: $(OUTPUT_DIR)/Calculator_Cosmology.o
$(OUTPUT_DIR)/CosmologyTypes.o: $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/likelihood.o $(OUTPUT_DIR)/GeneralTypes.o
$(OUTPUT_DIR)/Likelihood_Cosmology.o: $(OUTPUT_DIR)/Calculator_Cosmology.o $(OUTPUT_DIR)/CosmoTheory.o
$(OUTPUT_DIR)/CMBlikes.o: $(OUTPUT_DIR)/Likelihood_Cosmology.o
$(OUTPUT_DIR)/CMB.o: $(OUTPUT_DIR)/CMBlikes.o $(PLANCKLIKEFILES)
$(OUTPUT_DIR)/bbn.o: $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/likelihood.o $(OUTPUT_DIR)/Interpolation.o
$(OUTPUT_DIR)/mpk.o: $(OUTPUT_DIR)/CosmoTheory.o $(OUTPUT_DIR)/Likelihood_Cosmology.o
$(OUTPUT_DIR)/wigglez.o: $(OUTPUT_DIR)/mpk.o
$(OUTPUT_DIR)/bao.o: $(OUTPUT_DIR)/Likelihood_Cosmology.o
$(OUTPUT_DIR)/HST.o: $(OUTPUT_DIR)/Likelihood_Cosmology.o
$(OUTPUT_DIR)/cliklike.o: $(OUTPUT_DIR)/Likelihood_Cosmology.o
$(OUTPUT_DIR)/CosmologyConfig.o: $(PARAMETERIZATION) $(OUTPUT_DIR)/GeneralSetup.o $(OUTPUT_DIR)/Likelihood_Cosmology.o

$(OUTPUT_DIR)/supernovae_Union2.o: $(OUTPUT_DIR)/Likelihood_Cosmology.o
$(OUTPUT_DIR)/supernovae_SNLS.o: $(OUTPUT_DIR)/Likelihood_Cosmology.o
$(OUTPUT_DIR)/supernovae_JLA.o: $(OUTPUT_DIR)/Likelihood_Cosmology.o

$(OUTPUT_DIR)/Calculator_PICO.o: $(OUTPUT_DIR)/Calculator_CAMB.o $(OUTPUT_DIR)/fpico_interface.o


ifneq ($(highL),)
$(OUTPUT_DIR)/cliklike_CamSpec.o: $(highL)/libhigh.a
endif


ifneq ($(WMAP),)
$(OUTPUT_DIR)/CMB.o: $(WMAP)/libwmap9.a
F90FLAGS += $(IFLAG)$(cfitsio)/include $(IFLAG)$(WMAP) -DWMAP
LINKFLAGS +=  -L$(cfitsio)/lib -L$(WMAP) -lcfitsio -lwmap9
endif

F90FLAGS += $(MODOUT) $(IFLAG)$(OUTPUT_DIR)/

export FFLAGS
export F90C
export OUTPUT_DIR

directories:
	mkdir -p $(OUTPUT_DIR)

$(OUTPUT_DIR)/%.o: %.c
	$(CC) $(GSLINC) -c $*.c -o $(OUTPUT_DIR)/$*.o

$(OUTPUT_DIR)/%.o: %.f90 Makefile
	$(F90C) $(F90FLAGS) -c $*.f90 -o $(OUTPUT_DIR)/$*.o

$(OUTPUT_DIR)/%.o: %.F90 Makefile
	$(F90C) $(F90FLAGS) -c $*.F90 -o $(OUTPUT_DIR)/$*.o

$(OUTPUT_DIR)/fpico_interface.o: $(PICO)/fpico_interface.f90 Makefile
	$(F90C) $(F90FLAGS) -c $(PICO)/fpico_interface.f90 -o $(OUTPUT_DIR)/fpico_interface.o

cosmomc: directories camb $(OBJFILES)
	$(F90C) -o ../cosmomc $(OBJFILES) $(LINKFLAGS) $(F90FLAGS)

cosmomc_debug: directories camb $(OBJFILES)
	$(F90C) -o ../cosmomc_debug $(OBJFILES) $(LINKFLAGS) $(F90FLAGS)


clean: cleancosmomc
	rm -f ../camb/*.o ../camb/*.obj ../camb/*.mod ../camb/libcamb_$(RECOMBINATION).a

cleancosmomc:
	rm -f $(OUTPUT_DIR)/*.o $(OUTPUT_DIR)/*.mod ./Debug/* ./Release/*  ./DebugMPI/* ./ReleaseMPI/* ../core


getdist: directories $(DISTFILES)
	$(F90C) -o ../getdist $(DISTFILES) $(LAPACKL) $(GETDISTFLAGS)

#If you want to use the ppf module, change equations_ppf to equations
#and halofit_ppf to halofit also comment out wa_ppf in CMB_Cls_simple.f90
camb:
	cd ../camb && \
	$(MAKE) --file=Makefile_main libcamb_$(RECOMBINATION).a \
	RECOMBINATION=$(RECOMBINATION) EQUATIONS=equations_ppf NONLINEAR=halofit_ppf

$(highL)/libhigh.a:
	cd $(highL); make libhigh.a;

$(WMAP)/libwmap9.a:
	cd $(WMAP); make libwmap9.a;
