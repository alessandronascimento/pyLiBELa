#Updated version
#Name of the project
PROJ_NAME=pyLiBELa

#Classes

CLASSES= PARSER    \
	 Mol2      \
	 WRITER    \
	 Grid      \
	 COORD_MC  \
	 FindHB    \
	 Energy2   \
	 Conformer \
	 RAND      \
	 Gaussian  \
	 McEntropy \
	 SA	   \
	 MC	   \
	 Optimizer \
	 Docker    \
	 FullSearch

#Sources
SOURCE=$(wildcard src/*.cpp)

#Objects
OBJ=$(SOURCE:.cpp=.o)

#Compiler
CC=g++

#Flags for compiler
CC_FLAGS=-fPIC                       \
         -O3                         \
	 -I/usr/include/python3.10   \
	 -I/usr/include/openbabel3   \
	 -I/usr/include/eigen3       \
	 -DBUILD=0                   \
	 -DHAVE_EIGEN                \
	 -Wunused-value
	 
#Linkage
LD_FLAGS=-L/usr/lib/x86_64-linux-gnu \
         -lboost_python310            \
	 -lz                         \
	 -lopenbabel                 \
	 -lgsl                       \
	 -lnlopt_cxx                 \
	 -lgsl                       \
	 -lgslcblas                  \
	 -lm 


#Command to clean target
RM = rm -rf


#Compilation and linking

all: $(CLASSES)
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/rec.mol2.gz
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/lig.mol2.gz
	@ mkdir test
	@ mv *.gz test/
	@ echo ' '

FILES:
	@ git clone https://github.com/alessandronascimento/pyLiBELa.git
	@ mv pyLiBELa/src src
	@ rm -rf pyLiBELa
	@ mkdir -p obj

PARSER: OBJECTS 
	@echo 'Building PARSER dynamic library.'
	$(CC) -shared obj/pyPARSER.o -o pyPARSER.so $(LD_FLAGS)
	@echo ' '

Mol2:   OBJECTS
	@echo 'Building MOL2 dynamic library.'
	$(CC) -shared obj/pyMol2.o obj/pyPARSER.o -o pyMol2.so $(LD_FLAGS)
	@echo ' '
	
WRITER: OBJECTS
	@echo 'Building WRITER dynamic library.'
	$(CC) -shared obj/pyWRITER.o obj/pyMol2.o obj/pyPARSER.o -o pyWRITER.so $(LD_FLAGS)
	@echo ' '

Grid:   OBJECTS
	@echo 'Building GRID dynamic library.'
	$(CC) -shared obj/pyGrid.o obj/pyMol2.o obj/pyPARSER.o obj/pyWRITER.o -o pyGrid.so $(LD_FLAGS)
	@echo ' '

COORD_MC: OBJECTS
	@echo 'Building COORD_MC dynamic library.'
	$(CC) -shared obj/pyCOORD_MC.o obj/pyRAND.o obj/pyMol2.o -o pyCOORD_MC.so $(LD_FLAGS)
	@echo ' '
	
FindHB: OBJECTS
	@echo 'Building FIND_HB dynamic library.'
	$(CC) -shared obj/pyFindHB.o obj/pyMol2.o obj/pyPARSER.o obj/pyCOORD_MC.o -o pyFindHB.so $(LD_FLAGS)
	@echo ' '

Energy2: OBJECTS
	@echo 'Building ENERGY2 dynamic library.'
	$(CC) -shared obj/pyEnergy2.o obj/pyPARSER.o obj/pyMol2.o obj/pyGrid.o obj/pyWRITER.o -o pyEnergy2.so $(LD_FLAGS)
	@echo ' '

Conformer: OBJECTS
	@echo 'Building CONFORMER dynamic library.'
	$(CC) -shared obj/pyConformer.o obj/pyPARSER.o obj/pyMol2.o -o pyConformer.so $(LD_FLAGS)
	@echo ' '

RAND: OBJECTS
	@echo 'Building RAND dynamic library.'
	$(CC) -shared obj/pyRAND.o obj/pyPARSER.o obj/pyMol2.o -o pyRAND.so $(LD_FLAGS)
	@echo ' '

Gaussian: OBJECTS
	@echo 'Building GAUSSIAN dynamic library.'
	$(CC) -shared obj/pyGaussian.o obj/pyMol2.o obj/pyPARSER.o -o pyGaussian.so $(LD_FLAGS)
	@echo ' '

McEntropy: OBJECTS
	@echo 'Building MCENTROPY dynamic library.'
	$(CC) -shared obj/pyMcEntropy.o obj/pyPARSER.o obj/pyMol2.o obj/pyCOORD_MC.o -o pyMcEntropy.so $(LD_FLAGS)
	@echo ' '

SA: OBJECTS
	@echo 'Building SA dynamic library.'
	$(CC) -shared obj/pySA.o obj/pyMol2.o obj/pyGrid.o obj/pyEnergy2.o obj/pyCOORD_MC.o obj/pyWRITER.o obj/pyPARSER.o -o pySA.so $(LD_FLAGS)
	@echo ' '

MC: OBJECTS
	@echo 'Building MC dynamic library.'
	$(CC) -shared obj/pyMC.o obj/pyPARSER.o obj/pyGrid.o obj/pyMol2.o obj/pyCOORD_MC.o obj/pyEnergy2.o obj/pyWRITER.o obj/pyMcEntropy.o obj/pyOptimizer.o -o pyMC.so $(LD_FLAGS)
	@echo ' '

Optimizer: OBJECTS
	@echo 'Building OPTIMIZER dynamic library.'
	$(CC) -shared obj/pyOptimizer.o obj/pyMol2.o obj/pyEnergy2.o obj/pyPARSER.o obj/pyCOORD_MC.o obj/pyGaussian.o obj/pyWRITER.o obj/pySA.o obj/pyGrid.o -o pyOptimizer.so $(LD_FLAGS)
	@echo ' '

Docker: OBJECTS
	@echo 'Building DOCKER dynamic library.'
	$(CC) -shared obj/pyMol2.o obj/pyOptimizer.o obj/pyPARSER.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o -o pyDocker.so $(LD_FLAGS)
	@echo ' '

FullSearch: OBJECTS
	@echo 'Building FULLSEARCH dynamic library.'
	$(CC) -shared obj/pyMol2.o obj/pyPARSER.o obj/pyConformer.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o obj/pyEnergy2.o -o pyFullSearch.so $(LD_FLAGS)
	@echo ' '
	
OBJECTS: $(OBJ)
	@ echo 'Objects are created!'
	@ mv *.o obj

.cpp.o:
	$(CC) $(CC_FLAGS) $< -c
	@ echo ''

clean:
	@ $(RM) obj src test $(PROJ_NAME) *~
	@ $(RM) *.so *.o *.h *.cpp

.PHONY: all clean
