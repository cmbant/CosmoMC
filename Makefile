default: all

cosmomc: BUILD ?= MPI
cosmomc_debug: BUILD ?= MPI

getdist: ./source/*.*90
	cd ./source && make getdist BUILD=$(BUILD)

cosmomc: ./source/*.*90
	cd ./source && make cosmomc BUILD=$(BUILD)

cosmomc_debug: ./source/*.*90
	cd ./source && make cosmomc_debug OUTPUT_DIR=Debug BUILD=$(BUILD)

camspec: ./source/*.*90
	cd ./source && make highL=../highL PLANCKLIKE=cliklike_CamSpec

clean:
	cd ./source && make clean

all: cosmomc getdist
