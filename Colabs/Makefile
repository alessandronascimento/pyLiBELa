#Updated version
#Name of the project
PROJ_NAME=pyLiBELa

#Sources
SOURCE=$(wildcard src/*.cpp)

#Objects
OBJ=$(SOURCE:.cpp=.o)

#Compiler
CC=g++

#Flags for compiler
CC_FLAGS=-fPIC                       \
         -O3                         \
	 -I/usr/include/python3.10 \
	 -I/usr/local/include/openbabel3   \
	 -I/usr/include/eigen3       \
	 -DBUILD=0                   \
	 -DHAVE_EIGEN                \
	 -Wunused-value
	 
#Linkage
LD_FLAGS=-L/usr/lib/x86_64-linux-gnu \
         -lboost_python310           \
	 -lz                         \
	 -lopenbabel                 \
	 -lnlopt_cxx                 \
	 -lgsl                       \
	 -lgslcblas                  \
	 -lm 


NVFLAGS=-DBUILD=0 -DHAVE_EIGEN -I/usr/include/openbabel3 -I/usr/include/eigen3 -I/usr/include/python3.10

CUDAFLAGS=-L/usr/local/cuda/lib64 -lcudart

#Command to clean target
RM = rm -rf

#Compilation and linking

all: pyPARSER.so pyMol2.so pyWRITER.so pyGrid.so pyCOORD_MC.so pyFindHB.so pyEnergy2.so pyConformer.so pyRAND.so pyGaussian.so pyMcEntropy.so pySA.so pyMC.so pyOptimizer.so pyDocker.so pyFullSearch.so pyEngine.so McCharge smiles2mol2
#	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/rec.mol2.gz
#	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/lig.mol2.gz
#	@ mkdir test
#	@ mv *.gz test/
	@ echo ' '

cudaGrid: 
	@echo 'Building CUDA objects'
	nvcc --compiler-options '-fPIC' -c src/cudaGrid.cu -o obj/cudaGrid.o $(NVFLAGS)
	@echo ' '

pyPARSER.so: obj/pyPARSER.o
	@echo 'Building PARSER dynamic library.'
	$(CC) -shared obj/pyPARSER.o -o pyPARSER.so $(LD_FLAGS)
	@echo ' '

pyMol2.so: obj/pyMol2.o obj/pyPARSER.o
	@echo 'Building MOL2 dynamic library.'
	$(CC) -shared obj/pyMol2.o obj/pyPARSER.o -o pyMol2.so $(LD_FLAGS)
	@echo ' '
	
pyWRITER.so : obj/pyWRITER.o obj/pyMol2.o obj/pyPARSER.o
	@echo 'Building WRITER dynamic library.'
	$(CC) -shared obj/pyWRITER.o obj/pyMol2.o obj/pyPARSER.o -o pyWRITER.so $(LD_FLAGS)
	@echo ' '

pyGrid.so:   obj/pyGrid.o obj/pyMol2.o obj/pyPARSER.o obj/pyWRITER.o cudaGrid
	@echo 'Building GRID dynamic library.'
	$(CC) -shared obj/pyGrid.o obj/pyMol2.o obj/pyPARSER.o obj/pyWRITER.o obj/cudaGrid.o -o pyGrid.so $(LD_FLAGS) $(CUDAFLAGS)
	@echo ' '

pyCOORD_MC.so : obj/pyCOORD_MC.o obj/pyRAND.o obj/pyMol2.o
	@echo 'Building COORD_MC dynamic library.'
	$(CC) -shared obj/pyCOORD_MC.o obj/pyRAND.o obj/pyMol2.o -o pyCOORD_MC.so $(LD_FLAGS)
	@echo ' '
	
pyFindHB.so: obj/pyFindHB.o obj/pyMol2.o obj/pyPARSER.o obj/pyCOORD_MC.o
	@echo 'Building FIND_HB dynamic library.'
	$(CC) -shared obj/pyFindHB.o obj/pyMol2.o obj/pyPARSER.o obj/pyCOORD_MC.o -o pyFindHB.so $(LD_FLAGS)
	@echo ' '

pyEnergy2.so : obj/pyEnergy2.o obj/pyPARSER.o obj/pyMol2.o obj/pyGrid.o obj/pyWRITER.o
	@echo 'Building ENERGY2 dynamic library.'
	$(CC) -shared obj/pyEnergy2.o obj/pyPARSER.o obj/pyMol2.o obj/pyGrid.o obj/cudaGrid.o obj/pyWRITER.o -o pyEnergy2.so $(LD_FLAGS) $(CUDAFLAGS)
	@echo ' '

pyConformer.so: obj/pyConformer.o obj/pyPARSER.o obj/pyMol2.o
	@echo 'Building CONFORMER dynamic library.'
	$(CC) -shared obj/pyConformer.o obj/pyPARSER.o obj/pyMol2.o -o pyConformer.so $(LD_FLAGS)
	@echo ' '

pyRAND.so: obj/pyRAND.o obj/pyPARSER.o obj/pyMol2.o 
	@echo 'Building RAND dynamic library.'
	$(CC) -shared obj/pyRAND.o obj/pyPARSER.o obj/pyMol2.o -o pyRAND.so $(LD_FLAGS)
	@echo ' '

pyGaussian.so: obj/pyGaussian.o obj/pyMol2.o obj/pyPARSER.o 
	@echo 'Building GAUSSIAN dynamic library.'
	$(CC) -shared obj/pyGaussian.o obj/pyMol2.o obj/pyPARSER.o -o pyGaussian.so $(LD_FLAGS)
	@echo ' '

pyMcEntropy.so: obj/pyMcEntropy.o obj/pyPARSER.o obj/pyMol2.o obj/pyCOORD_MC.o
	@echo 'Building MCENTROPY dynamic library.'
	$(CC) -shared obj/pyMcEntropy.o obj/pyPARSER.o obj/pyMol2.o obj/pyCOORD_MC.o -o pyMcEntropy.so $(LD_FLAGS)
	@echo ' '

pySA.so: obj/pySA.o obj/pyMol2.o obj/pyGrid.o obj/pyEnergy2.o obj/pyCOORD_MC.o obj/pyWRITER.o obj/pyPARSER.o
	@echo 'Building SA dynamic library.'
	$(CC) -shared obj/pySA.o obj/pyMol2.o obj/pyGrid.o obj/pyEnergy2.o obj/pyCOORD_MC.o obj/pyWRITER.o obj/pyPARSER.o obj/cudaGrid.o -o pySA.so $(LD_FLAGS) $(CUDAFLAGS)
	@echo ' '

pyMC.so: obj/pyMC.o obj/pyPARSER.o obj/pyGrid.o obj/pyMol2.o obj/pyCOORD_MC.o obj/pyEnergy2.o obj/pyWRITER.o obj/pyMcEntropy.o obj/pyOptimizer.o obj/pyGaussian.o obj/pySA.o
	@echo 'Building MC dynamic library.'
	$(CC) -shared obj/pyMC.o obj/pyPARSER.o obj/pyGrid.o obj/pyMol2.o obj/pyCOORD_MC.o obj/pyEnergy2.o obj/pyWRITER.o obj/pyMcEntropy.o obj/pyOptimizer.o obj/pyGaussian.o obj/pySA.o obj/cudaGrid.o -o pyMC.so $(LD_FLAGS) $(CUDAFLAGS)
	@echo ' '

pyOptimizer.so: obj/pyOptimizer.o obj/pyMol2.o obj/pyEnergy2.o obj/pyPARSER.o obj/pyCOORD_MC.o obj/pyGaussian.o obj/pyWRITER.o obj/pySA.o obj/pyMC.o obj/pyGrid.o obj/pyMcEntropy.o
	@echo 'Building OPTIMIZER dynamic library.'
	$(CC) -shared obj/pyOptimizer.o obj/pyMol2.o obj/pyEnergy2.o obj/pyPARSER.o obj/pyCOORD_MC.o obj/pyGaussian.o obj/pyWRITER.o obj/pySA.o obj/pyMC.o obj/pyGrid.o obj/pyMcEntropy.o obj/cudaGrid.o -o pyOptimizer.so $(LD_FLAGS) $(CUDAFLAGS)
	@echo ' '

pyDocker.so: obj/pyDocker.o obj/pyMol2.o obj/pyOptimizer.o obj/pyPARSER.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o obj/pySA.o obj/pyEnergy2.o obj/pyGaussian.o 
	@echo 'Building DOCKER dynamic library.'
	$(CC) -shared obj/pyDocker.o obj/pyMol2.o obj/pyOptimizer.o obj/pyPARSER.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o obj/pySA.o obj/pyEnergy2.o obj/pyGaussian.o obj/cudaGrid.o -o pyDocker.so $(LD_FLAGS) $(CUDAFLAGS)
	@echo ' '

pyFullSearch.so: obj/pyFullSearch.o obj/pyMol2.o obj/pyPARSER.o obj/pyConformer.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o obj/pyEnergy2.o
	@echo 'Building FULLSEARCH dynamic library.'
	$(CC) -shared obj/pyFullSearch.o obj/pyMol2.o obj/pyPARSER.o obj/pyConformer.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o obj/pyEnergy2.o obj/cudaGrid.o -o pyFullSearch.so $(LD_FLAGS) $(CUDAFLAGS)
	@echo ' '

pyEngine.so: obj/pyEngine.o obj/pyMol2.o obj/pyPARSER.o obj/pyConformer.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o obj/pyDocker.o obj/pyOptimizer.o obj/pySA.o obj/pyEnergy2.o obj/pyGaussian.o obj/pyFindHB.o
	@echo 'Building pyEngine dynamic library.'
	$(CC) -shared obj/pyEngine.o obj/pyMol2.o obj/pyPARSER.o obj/pyConformer.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o obj/pyDocker.o obj/pyOptimizer.o obj/pySA.o obj/pyEnergy2.o obj/pyGaussian.o obj/pyFindHB.o -o pyEngine.so $(LD_FLAGS)
	@echo ' '

McCharge: obj/pyMol2.o obj/pyPARSER.o
	@echo 'Building McCharge tool...'
	$(CC) $(CC_FLAGS) test/McCharge.cpp -c -o test/McCharge.o 
	$(CC) obj/pyMol2.o obj/pyPARSER.o test/McCharge.o -o test/McCharge $(LD_FLAGS) -lpython3.10
	@echo ' '

smiles2mol2: obj/pyMol2.o
	@echo 'Building smiles2mol2 tool...'
	$(CC) $(CC_FLAGS) test/smiles2mol2.cpp -c -o test/smiles2mol2.o 
	$(CC) obj/pyMol2.o test/smiles2mol2.o -o test/smiles2mol2 $(LD_FLAGS) -lpython3.10
	@echo ' '
	
obj/%.o: src/%.cpp
	$(CC) $(CC_FLAGS) $< -c -o $@
	@ echo ''

clean:
	@ $(RM) obj src test $(PROJ_NAME) *~
	@ $(RM) *.so *.o *.h *.cpp

.PHONY: all clean
