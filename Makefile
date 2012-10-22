default: all

getdist: ./source/*.*90
	cd ./source && make getdist

cosmomc: ./source/*.*90
	cd ./source && make

all: ./source/*.*90
	cd ./source && make all

