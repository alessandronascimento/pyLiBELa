#!/bin/bash

CC=g++
CFLAGS=-c -O3 -ffast-math -fopenmp -DBUILD=0 -DHAVE_EIGEN -I/usr/include/openbabel3 -I/usr/include/eigen3 -I/usr/include/python3.11
NVFLAGS=-DBUILD=0 -DHAVE_EIGEN -I/usr/include/openbabel3 -I/usr/include/eigen3 -I/usr/include/python3.11
LDFLAGS=-lz -L/usr/lib/x86_64-linux-gnu -lopenbabel -lgomp -lpthread -lgsl -lgslcblas -lm -lpython3.11 -lboost_python311 -L/usr/local/cuda/lib64 -lcudart
SOURCES=McGrid.cpp 
OBJECTS=$(SOURCES:.cpp=.o) cudaGrid.o cudaEnergy2.o

all: McGrid

McConf: McConf.o
	$(CC) McConf.o -o McConf $(LDFLAGS)

cudaGrid:
	nvcc -c ../cudaGrid.o ../cudaGrid.cu $(NVFLAGS)

cudaEnergy2: 
	nvcc -c -arch=sm_72 ../cudaEnergy2.o ../cudaEnergy2.cu $(NVFLAGS)

McGrid: McGrid.o cudaGrid cudaEnergy2
	$(CC) McGrid.o -o McGrid cudaGrid.o cudaEnergy2.o  $(LDFLAGS) 

mol2_rmsd: mol2_rmsd.o
	$(CC) mol2_rmsd.o -o mol2_rmsd $(LDFLAGS)

MrPocket: MrPocket.o
	$(CC) MrPocket.o -o McPocket $(LDFLAGS)

McColony: McColony.o 
	$(CC) McColony.o -o McColony $(LDFLAGS)

McSummation: FESummation.o
	$(CC) FESummation.o -o McSummation $(LDFLAGS)

McBootstrapper: McBootstrapper.o
	$(CC) McBootstrapper.o -o McBootstrapper $(LDFLAGS)

McCharge: McCharge.o
	$(CC) McCharge.o -o McCharge $(LDFLAGS)

AutoCorr: AutoCorrelation.o
	$(CC) AutoCorrelation.o -o AutoCorr $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS) McConf McGrid mol2_rmsd McPocket
