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

#Sources
SOURCE=$(wildcard ./src/*.cpp)

#Objects
OBJ=$(subst .cpp,.o, $(subst src,obj,$(SOURCE)))

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

all: FILES $(PROJ_NAME)

$(PROJ_NAME): $(CLASSES)
	
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
	@ mv *.h *.cpp src/


./obj/%.o: ./src/%.cpp
	$(CC) -o $@ $< $(CC_FLAGS)


PARSER: $(OBJ)
	@echo 'Building PARSER dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyPARSER.o -o pyPARSER.so $(LD_FLAGS)
	@echo ' '

Mol2: pyMol2.o PARSER
	@echo 'Building MOL2 dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyMol2.o obj/pyPARSER.o -o pyMol2.so $(LD_FLAGS)
	@echo ' '

WRITER: Mol2 PARSER
	@echo 'Building WRITER dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyWRITER.o obj/pyMol2.o obj/pyPARSER.o -o pyWRITER.so $(LD_FLAGS)
	@echo ' '

Grid:   Mol2 PARSER WRITER
	@echo 'Building GRID dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyGrid.o -o pyGrid.so $(LD_FLAGS)
	@echo ' '

COORD_MC: RAND Mol2
	@echo 'Building COORD_MC dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyCOORD_MC.o -o pyCOORD_MC.so $(LD_FLAGS)
	@echo ' '
	
FindHB: Mol2 PARSER COORD_MC
	@echo 'Building FIND_HB dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyFindHB.o -o pyFindHB.so $(LD_FLAGS)
	@echo ' '

Energy2: PARSER Mol2 Grid
	@echo 'Building ENERGY2 dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyEnergy2.o -o pyEnergy2.so $(LD_FLAGS)
	@echo ' '

Conformer: PARSER Mol2
	@echo 'Building CONFORMER dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyConformer.o -o pyConformer.so $(LD_FLAGS)
	@echo ' '

RAND: PARSER Mol2
	@echo 'Building RAND dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyRAND.o -o pyRAND.so $(LD_FLAGS)
	@echo ' '

Gaussian: Mol2 PARSER
	@echo 'Building GAUSSIAN dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyGaussian.o -o pyGaussian.so $(LD_FLAGS)
	@echo ' '

McEntropy: PARSER Mol2 COORD_MC
	@echo 'Building MCENTROPY dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyMcEntropy.o -o pyMcEntropy.so $(LD_FLAGS)
	@echo ' '

SA: Mol2 Grid Energy2 COORD_MC PARSER
	@echo 'Building SA dynamic library.'
	$(CC) $(CC_FLAGS) obj/pySA.o -o pySA.so $(LD_FLAGS)
	@echo ' '

MC: PARSER Grid Mol2 COORD_MC Energy2 WRITER McEntropy 
	@echo 'Building MC dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyMC.o -o pyMC.so $(LD_FLAGS)
	@echo ' '

Optimizer: Mol2 Energy2 PARSER COORD_MC Gaussian WRITER SA MC
	@echo 'Building OPTIMIZER dynamic library.'
	$(CC) $(CC_FLAGS) obj/pyOptimizer.o -o pyOptimizer.so $(LD_FLAGS)
	@echo ' '

clean:
	@ $(RM) *.o src/*.h src/*.cpp $(PROJ_NAME) *~

.PHONY: all clean
