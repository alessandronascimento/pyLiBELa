#Name of the project
PROJ_NAME=pyLiBELa

#Classes
CLASS= PARSER Mol2 WRITER Grid COORD_MC FindHB Energy2 Conformer RAND
#Compiler
CC=g++

#Flags for compiler
CC_FLAGS=-fPIC                       \
	 -shared                     \
	 -I/usr/include/python3.9   \
	 -I/usr/include/openbabel3   \
	 -I/usr/include/eigen3       \
	 -DBUILD=0                   \
	 -DHAVE_EIGEN
	 
#Linkage
LD_FLAGS=-L/usr/lib/x86_64-linux-gnu \
	 -lboost_python38            \
	 -lz                         \
	 -lopenbabel                 \
	 -lgsl


#Command to clean target
RM = rm -rf


#Compilation and linking

all: $(PROJ_NAME)

$(PROJ_NAME): $(CLASSES)
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/rec.mol2.gz
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/lig.mol2.gz
#	@ mv *.h *.cpp src/
#	@ mv *gz test/

PARSER: 
	@echo 'Building PARSER dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyPARSER.cpp
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyPARSER.h
	$(CC) $(CC_FLAGS) pyPARSER.cpp -o pyPARSER.so $(LD_FLAGS)
	@echo ' '

Mol2:
	@echo 'Building Mol2 dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMol2.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMol2.cpp
	$(CC) $(CC_FLAGS) pyMol2.cpp -o pyMol2.so $(LD_FLAGS)
	@echo ' '

WRITER:
	@echo 'Building WRITER dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/iMcLiBELa.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyWRITER.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyWRITER.cpp
	$(CC) $(CC_FLAGS) pyWRITER.cpp -o pyWRITER.so $(LD_FLAGS)
	@echo ' '

Grid:
	@echo 'Building Grid dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyGrid.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyGrid.cpp
	$(CC) $(CC_FLAGS) pyGrid.cpp -o pyGrid.so $(LD_FLAGS)
	@echo ' '

COORD_MC:
	@echo 'Building COORD_MC dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyCOORD_MC.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyCOORD_MC.cpp
	$(CC) $(CC_FLAGS) pyCOORD_MC.cpp -o pyCOORD_MC.so $(LD_FLAGS)
	@echo ' '

FindHB:
	@echo 'Building FindHB dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyFindHB.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyFindHB.cpp
	$(CC) $(CC_FLAGS) pyFindHB.cpp -o pyFindHB.so $(LD_FLAGS)
	@echo ' '

Energy2:
	@echo 'Building Energy2 dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyEnergy2.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyEnergy2.cpp
	$(CC) $(CC_FLAGS) pyEnergy2.cpp -o pyEnergy2.so $(LD_FLAGS)
	@echo ' '

Conformer:
	@echo 'Building Conformer dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyConformer.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyConformer.cpp
	$(CC) $(CC_FLAGS) pyConformer.cpp -o pyConformer.so $(LD_FLAGS)
	@echo ' '

RAND:
	@echo 'Building RAND dinamic library.'
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyRAND.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyRAND.cpp
	$(CC) $(CC_FLAGS) pyRAND.cpp -o pyRAND.so $(LD_FLAGS)
	@echo ' '

clean:
	@ $(RM) *.so *.h *.cpp $(PROJ_NAME) *~

.PHONY: all clean
