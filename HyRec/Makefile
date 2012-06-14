CC = icc
AR = ar rvs
CCFLAG = -O2 -DCAMB -DHYRECPATH=\"$(PWD)/\" 
LDFLAG = -O2

%.o: %.c 
	$(CC) $(CCFLAG) -c $*.c -o $*.o

HYREC_SRC = hyrectools.o helium.o hydrogen.o history.o 
HYREC_EXE = hyrec.o 

default: libhyrec.a

clean: 
	rm *.o

hyrec: $(HYREC_SRC) $(HYREC_EXE)
	$(CC) $(LDFLAG) -o hyrec $(HYREC_SRC) $(HYREC_EXE) -lm

libhyrec.a: $(HYREC_SRC)
	$(AR) $@ $(HYREC_SRC)
