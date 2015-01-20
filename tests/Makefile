OUTPUT_DIR=../source/Debug

UTILS = $(OUTPUT_DIR)/MiscUtils.o $(OUTPUT_DIR)/ArrayUtils.o $(OUTPUT_DIR)/StringUtils.o $(OUTPUT_DIR)/MpiUtils.o $(OUTPUT_DIR)/FileUtils.o \
 $(OUTPUT_DIR)/RandUtils.o $(OUTPUT_DIR)/ObjectLists.o $(OUTPUT_DIR)/Interpolation.o

default:
	gfortran -J$(OUTPUT_DIR) $(UTILS) ObjectLists_tests.f90 interpolation_tests.f90 general_unittests.f90