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
	@ mkdir src
	@ mv *.gz test/
	@ mv *.h *.cpp src/

SOURCES:
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

PARSER: 
	@echo 'Building PARSER dinamic library.'
	$(CC) $(CC_FLAGS) pyPARSER.cpp -o pyPARSER.so $(LD_FLAGS)
	@echo ' '

Mol2: PARSER
	@echo 'Building Mol2 dinamic library.'
	$(CC) $(CC_FLAGS) pyMol2.cpp -o pyMol2.so $(LD_FLAGS)
	@echo ' '

WRITER: Mol2 PARSER
	@echo 'Building WRITER dinamic library.'
	$(CC) $(CC_FLAGS) pyWRITER.cpp -o pyWRITER.so $(LD_FLAGS)
	@echo ' '

Grid:   Mol2 PARSER WRITER
	@echo 'Building Grid dinamic library.'
	$(CC) $(CC_FLAGS) pyGrid.cpp -o pyGrid.so $(LD_FLAGS)
	@echo ' '

COORD_MC: RAND Mol2
	@echo 'Building COORD_MC dinamic library.'
	$(CC) $(CC_FLAGS) pyCOORD_MC.cpp -o pyCOORD_MC.so $(LD_FLAGS)
	@echo ' '

FindHB: Mol2 PARSER COORD_MC
	@echo 'Building FindHB dinamic library.'
	$(CC) $(CC_FLAGS) pyFindHB.cpp -o pyFindHB.so $(LD_FLAGS)
	@echo ' '

Energy2: PARSER Mol2 Grid
	@echo 'Building Energy2 dinamic library.'
	$(CC) $(CC_FLAGS) pyEnergy2.cpp -o pyEnergy2.so $(LD_FLAGS)
	@echo ' '

Conformer: PARSER Mol2
	@echo 'Building Conformer dinamic library.'
	$(CC) $(CC_FLAGS) pyConformer.cpp -o pyConformer.so $(LD_FLAGS)
	@echo ' '

RAND: PARSER Mol2
	@echo 'Building RAND dinamic library.'
	$(CC) $(CC_FLAGS) pyRAND.cpp -o pyRAND.so $(LD_FLAGS)
	@echo ' '

Gaussian: Mol2 PARSER
	@echo 'Building Gaussian dinamic library.'
	$(CC) $(CC_FLAGS) pyGaussian.cpp -o pyGaussian.so $(LD_FLAGS)
	@echo ' '

McEntropy: PARSER Mol2 COORD_MC
	@echo 'Building McEntropy dinamic library.'
	$(CC) $(CC_FLAGS) pyMcEntropy.cpp -o pyMcEntropy.so $(LD_FLAGS)
	@echo ' '

SA: Mol2 Grid Energy2 COORD_MC PARSER
	@echo 'Building SA dinamic library.'
	$(CC) $(CC_FLAGS) pySA.cpp -o pySA.so $(LD_FLAGS)
	@echo ' '

MC: PARSER Grid Mol2 COORD_MC Energy2 WRITER McEntropy 
	@echo 'Building MC dinamic library.'
	$(CC) $(CC_FLAGS) pyMC.cpp -o pyMC.so $(LD_FLAGS)

	@echo ' '

Optimizer: Mol2 Energy2 PARSER COORD_MC Gaussian WRITER SA MC
	$(CC) $(CC_FLAGS) pyOptimizer.cpp -o pyOptimizer.so $(LD_FLAGS)

	@echo ' '


clean:
	@ $(RM) *.so src/*.h src/*.cpp $(PROJ_NAME) *~

.PHONY: all clean
