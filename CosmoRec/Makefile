#============================================================================================
# definitions
#============================================================================================
DEV_DIR = ./Development
REC_BASE= ./Rec_database

MAKE_COMMAND = make -f Makefile.in

#============================================================================================
# rules
#============================================================================================
all: 
	$(MAKE_COMMAND) all

lib: 
	$(MAKE_COMMAND) lib

bin: 
	$(MAKE_COMMAND) bin

#============================================================================================
# rules to clean up
#============================================================================================
clean:
	rm -f *.o
	rm -f ./Modules/*.o 

cleanPDE:
	rm -f ./PDE_Problem/*.o 

cleanall:
	rm -f *.o *~ CosmoRec libCosmoRec.a
	rm -f ./Modules/*.o ./Modules/*~
	rm -f ./PDE_Problem/*.o ./PDE_Problem/*~ 
	rm -f $(REC_BASE)/Effective_Rates.HI/*.o $(REC_BASE)/Effective_Rates.HI/*~ 
	rm -f $(REC_BASE)/Effective_Rates.HeI/*.o $(REC_BASE)/Effective_Rates.HeI/*~ 

cleanallDEV:
	make cleanall
	rm -f $(DEV_DIR)/*.o $(DEV_DIR)/*.~
	rm -f $(DEV_DIR)/Definitions/*.o $(DEV_DIR)/Definitions/*.~
	rm -f $(DEV_DIR)/Simple_routines/*.o $(DEV_DIR)/Simple_routines/*.~
	rm -f $(DEV_DIR)/Cosmology/*.o $(DEV_DIR)/Cosmology/*.~
	rm -f $(DEV_DIR)/Recombination/*.o $(DEV_DIR)/Recombination/*~
	rm -f $(DEV_DIR)/Hydrogenic/*.o $(DEV_DIR)/Hydrogenic/*~
	rm -f $(DEV_DIR)/Line_profiles/*.o $(DEV_DIR)/Line_profiles/*~
	rm -f $(DEV_DIR)/Helium/*.o $(DEV_DIR)/Helium/*~
	rm -f $(DEV_DIR)/Integration/*.o $(DEV_DIR)/Integration/*~
	rm -f $(DEV_DIR)/ODE_PDE_Solver/*.o $(DEV_DIR)/ODE_PDE_Solver/*~
	cd $(DEV_DIR)/Cosmology/Recfast++/; make tidy;

tidy:
	make cleanallDEV

wipeDS:
	find . -type f -name \.DS_Store -print | xargs rm
