#Makefile --- High ell Likelihood code...

HIGHELLLIB = libhigh.a
OBJS =  Highell_options.o \
	Highell_subroutines.o \
	Foregrounds_loading.o \
	ACT_south_likelihood.o \
	ACT_equa_likelihood.o \
	SPT_keisler_likelihood.o \
	SPT_reichardt_likelihood.o \
	Highell_likelihood.o

RM = rm -f

LIBS ?= $(LAPACKL)


F90C    ?= ifort
FFLAGS  ?=  -O2 -ip -W0 -WB -fpp2 -vec_report0
INCS    ?= -I. -I/usr/local/include
LIBS  ?= -L. -L/usr/local/shared/lib -lmkl -lmkl_core -lmkl_lapack -lguide -lpthread -lm -llapack

#			Rules.
#
all: $(HIGHELLLIB) test

check: test
	./test

$(HIGHELLLIB): $(OBJS)
	ar r $@ $^

test: $(OBJS) test.o
	$(F90C) $(FFLAGS) -o $@ test.o $(OBJS) $(LIBS)

%.o: %.f90
	$(F90C) $(FFLAGS) $(INCS) -c -o $@ $<

%.o: %.F90
	$(F90C) $(FFLAGS) $(INCS) -c -o $@ $<

clean:
	$(RM) *.o *.mod *.log *~ *.a test

distclean: clean
	$(RM) test

