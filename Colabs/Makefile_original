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
	 Optimizer

#Compiler
CC=g++

#Flags for compiler
CC_FLAGS=-fPIC                       \
         -O3                         \
	 -shared                     \
	 -I/usr/include/python3.9    \
	 -I/usr/include/openbabel3   \
	 -I/usr/include/eigen3       \
	 -DBUILD=0                   \
	 -DHAVE_EIGEN
	 
#Linkage
LD_FLAGS=-L/usr/lib/x86_64-linux-gnu \
	 -lboost_python38            \
	 -lz                         \
	 -lopenbabel                 \
	 -lgsl                       \
	 -lnlopt_cxx


#Command to clean target
RM = rm -rf


#Compilation and linking

all:  $(PROJ_NAME)

$(PROJ_NAME): SOURCES $(CLASSES)
	
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/rec.mol2.gz
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/lig.mol2.gz
	@ mkdir test
	@ mv *.gz test/


SOURCES:
	@ mkdir src ;
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyPARSER.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyPARSER.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMol2.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMol2.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyWRITER.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyWRITER.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyGrid.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyGrid.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyCOORD_MC.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyCOORD_MC.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyFindHB.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyFindHB.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyEnergy2.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyEnergy2.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyConformer.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyConformer.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyRAND.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyRAND.cpp	
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyGaussian.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyGaussian.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMcEntropy.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMcEntropy.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/iMcLiBELa.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pySA.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pySA.cpp	
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMC.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMC.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyOptimizer.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyOptimizer.cpp
	@ mv *.h *.cpp src/

PARSER: 
	@echo 'Building PARSER dynamic library.'
	$(CC) $(CC_FLAGS) src/pyPARSER.cpp -o pyPARSER.so $(LD_FLAGS)
	@echo ' '

Mol2: PARSER
	@echo 'Building MOL2 dynamic library.'
	$(CC) $(CC_FLAGS) src/pyMol2.cpp -o pyMol2.so $(LD_FLAGS)
	@echo ' '

WRITER: Mol2 PARSER
	@echo 'Building WRITER dynamic library.'
	$(CC) $(CC_FLAGS) src/pyWRITER.cpp -o pyWRITER.so $(LD_FLAGS)
	@echo ' '

Grid:   Mol2 PARSER WRITER
	@echo 'Building GRID dynamic library.'
	$(CC) $(CC_FLAGS) src/pyGrid.cpp -o pyGrid.so $(LD_FLAGS)
	@echo ' '

COORD_MC: RAND Mol2
	@echo 'Building COORD_MC dynamic library.'
	$(CC) $(CC_FLAGS) src/pyCOORD_MC.cpp -o pyCOORD_MC.so $(LD_FLAGS)
	@echo ' '

FindHB: Mol2 PARSER COORD_MC
	@echo 'Building FIND_HB dynamic library.'
	$(CC) $(CC_FLAGS) src/pyFindHB.cpp -o pyFindHB.so $(LD_FLAGS)
	@echo ' '

Energy2: PARSER Mol2 Grid
	@echo 'Building ENERGY2 dynamic library.'
	$(CC) $(CC_FLAGS) src/pyEnergy2.cpp -o pyEnergy2.so $(LD_FLAGS)
	@echo ' '

Conformer: PARSER Mol2
	@echo 'Building CONFORMER dynamic library.'
	$(CC) $(CC_FLAGS) src/pyConformer.cpp -o pyConformer.so $(LD_FLAGS)
	@echo ' '

RAND: PARSER Mol2
	@echo 'Building RAND dynamic library.'
	$(CC) $(CC_FLAGS) src/pyRAND.cpp -o pyRAND.so $(LD_FLAGS)
	@echo ' '

Gaussian: Mol2 PARSER
	@echo 'Building GAUSSIAN dynamic library.'
	$(CC) $(CC_FLAGS) src/pyGaussian.cpp -o pyGaussian.so $(LD_FLAGS)
	@echo ' '

McEntropy: PARSER Mol2 COORD_MC
	@echo 'Building MCENTROPY dynamic library.'
	$(CC) $(CC_FLAGS) src/pyMcEntropy.cpp -o pyMcEntropy.so $(LD_FLAGS)
	@echo ' '

SA: Mol2 Grid Energy2 COORD_MC PARSER
	@echo 'Building SA dynamic library.'
	$(CC) $(CC_FLAGS) src/pySA.cpp -o pySA.so $(LD_FLAGS)
	@echo ' '

MC: PARSER Grid Mol2 COORD_MC Energy2 WRITER McEntropy 
	@echo 'Building MC dynamic library.'
	$(CC) $(CC_FLAGS) src/pyMC.cpp -o pyMC.so $(LD_FLAGS)
	@echo ' '

Optimizer: Mol2 Energy2 PARSER COORD_MC Gaussian WRITER SA MC
	@echo 'Building OPTIMIZER dynamic library.'
	$(CC) $(CC_FLAGS) src/pyOptimizer.cpp -o pyOptimizer.so $(LD_FLAGS)
	@echo ' '

clean:
	@ $(RM) *.so src/*.h src/*.cpp $(PROJ_NAME) *~

.PHONY: all clean
