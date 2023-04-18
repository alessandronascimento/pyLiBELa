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
SOURCE=$(wildcard ./src/*.cpp)

#Objects
OBJ=$(subst .cpp,.o, $(subst src,.,$(SOURCE)))

#Compiler
CC=g++

#Flags for compiler
CC_FLAGS=-fPIC                       \
         -O3                         \
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
	 -lnlopt_cxx                 \
	 -lgsl                       \
	 -lgslcblas                  \
	 -lm 


#Command to clean target
RM = rm -rf


#Compilation and linking

all: FILES OBJECTS LINKAGE
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/rec.mol2.gz
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/lig.mol2.gz
	@ mkdir test
	@ mv *.gz test/
	@ echo ' '

FILES:
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
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyDocker.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyDocker.cpp
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyFullSearch.h
	@ wget -q https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyFullSearch.cpp
	@ echo 'Downloaded files'
	@ mkdir src 
	@ mv *.h *.cpp src/

PARSER: 
	@echo 'Building PARSER dynamic library.'
	$(CC) -shared obj/pyPARSER.o -o pyPARSER.so $(LD_FLAGS)
	@echo ' '

Mol2:   PARSER
	@echo 'Building MOL2 dynamic library.'
	$(CC) -shared obj/pyMol2.o obj/pyPARSER.o -o pyMol2.so $(LD_FLAGS)
	@echo ' '
	
WRITER: Mol2 PARSER
	@echo 'Building WRITER dynamic library.'
	$(CC) -shared obj/pyWRITER.o obj/pyMol2.o obj/pyPARSER.o -o pyWRITER.so $(LD_FLAGS)
	@echo ' '

Grid:   Mol2 PARSER WRITER
	@echo 'Building GRID dynamic library.'
	$(CC) -shared obj/pyGrid.o obj/pyMol2.o obj/pyPARSER.o obj/pyWRITER.o -o pyGrid.so $(LD_FLAGS)
	@echo ' '

COORD_MC: RAND Mol2
	@echo 'Building COORD_MC dynamic library.'
	$(CC) -shared obj/pyCOORD_MC.o obj/pyRAND.o obj/pyMol2.o -o pyCOORD_MC.so $(LD_FLAGS)
	@echo ' '
	
FindHB: Mol2 PARSER COORD_MC
	@echo 'Building FIND_HB dynamic library.'
	$(CC) -shared obj/pyFindHB.o obj/pyMol2.o obj/pyPARSER.o obj/pyCOORD_MC.o -o pyFindHB.so $(LD_FLAGS)
	@echo ' '

Energy2: PARSER Mol2 Grid
	@echo 'Building ENERGY2 dynamic library.'
	$(CC) -shared obj/pyEnergy2.o obj/pyPARSER.o obj/pyMol2.o obj/pyGrid.o obj/pyWRITER.o -o pyEnergy2.so $(LD_FLAGS)
	@echo ' '

Conformer: PARSER Mol2
	@echo 'Building CONFORMER dynamic library.'
	$(CC) -shared obj/pyConformer.o obj/pyPARSER.o obj/pyMol2.o -o pyConformer.so $(LD_FLAGS)
	@echo ' '

RAND: PARSER Mol2
	@echo 'Building RAND dynamic library.'
	$(CC) -shared obj/pyRAND.o obj/pyPARSER.o obj/pyMol2.o -o pyRAND.so $(LD_FLAGS)
	@echo ' '

Gaussian: Mol2 PARSER
	@echo 'Building GAUSSIAN dynamic library.'
	$(CC) -shared obj/pyGaussian.o obj/pyMol2.o obj/pyPARSER.o -o pyGaussian.so $(LD_FLAGS)
	@echo ' '

McEntropy: PARSER Mol2 COORD_MC
	@echo 'Building MCENTROPY dynamic library.'
	$(CC) -shared obj/pyMcEntropy.o obj/pyPARSER.o obj/pyMol2.o obj/pyCOORD_MC.o -o pyMcEntropy.so $(LD_FLAGS)
	@echo ' '

SA: Mol2 Grid Energy2 COORD_MC PARSER
	@echo 'Building SA dynamic library.'
	$(CC) -shared obj/pySA.o obj/pyMol2.o obj/pyGrid.o obj/pyEnergy2.o obj/pyCOORD_MC.o obj/pyWRITER.o obj/pyPARSER.o -o pySA.so $(LD_FLAGS)
	@echo ' '

MC: PARSER Grid Mol2 COORD_MC Energy2 WRITER McEntropy #Optimizer 
	@echo 'Building MC dynamic library.'
	$(CC) -shared obj/pyMC.o obj/pyPARSER.o obj/pyGrid.o obj/pyMol2.o obj/pyCOORD_MC.o obj/pyEnergy2.o obj/pyWRITER.o obj/pyMcEntropy.o obj/pyOptimizer.o -o pyMC.so $(LD_FLAGS)
	@echo ' '

Optimizer: Mol2 Energy2 PARSER COORD_MC Gaussian WRITER SA MC
	@echo 'Building OPTIMIZER dynamic library.'
	$(CC) -shared obj/pyOptimizer.o obj/pyMol2.o obj/pyEnergy2.o obj/pyPARSER.o obj/pyCOORD_MC.o obj/pyGaussian.o obj/pyWRITER.o obj/pySA.o obj/pyMC.o -o pyOptimizer.so $(LD_FLAGS)
	@echo ' '

Docker: Mol2 Optimizer PARSER WRITER COORD_MC Grid
	@echo 'Building DOCKER dynamic library.'
	$(CC) -shared obj/pyMol2.o obj/pyOptimizer.o obj/pyPARSER.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o -o pyDocker.so $(LD_FLAGS)
	@echo ' '

FullSearch: Mol2 PARSER Conformer WRITER COORD_MC Grid Energy2
	@echo 'Building FULLSEARCH dynamic library.'
	$(CC) -shared obj/pyMol2.o obj/pyPARSER.o obj/pyConformer.o obj/pyWRITER.o obj/pyCOORD_MC.o obj/pyGrid.o obj/pyEnergy2.o -o pyFullSearch.so $(LD_FLAGS)
	@echo ' '
	
OBJECTS: $(OBJ)
	@ mkdir obj
	@ mv *.o obj/
	@ echo 'Succesfully build .o files!'

LINKAGE: $(CLASSES)
	@echo 'Classes are linked!'
	@echo ''

%.o: ./src/%.cpp
	@$(CC) $< $(CC_FLAGS) -c

clean:
	@ $(RM) obj src $(PROJ_NAME) *~

.PHONY: all clean
