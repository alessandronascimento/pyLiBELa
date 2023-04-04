#Name of the project
PROJ_NAME=pyLiBELa

#Every .cpp file
CPP_SRC= pyPARSER.cpp pyMol2.cpp


#Every .h file
H_SRC=$(wildcard ./src/*.h)

#Dinamic libraries
OBJ=$(subst .cpp,.so,$(subst src,obj,$(CPP_SRC)))


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
	 -lboost_python38           \
	 -lz                         \
	 -lopenbabel                  


#Command to clean target
RM = rm -rf


#Compilation and linking

all: objFolder $(PROJ_NAME)

$(PROJ_NAME): $(pyPARSER.so pyMol2.so)


pyPARSER.so: pyPARSER.cpp pyPARSER.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyPARSER.cpp
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyPARSER.h
	@echo 'Building PARSER dinamic library.'
	$(CC) $< $(CC_FLAGS) pyPARSER.cpp -o pyPARSER.so $(LD_FLAGS)
	@echo ' '

pyMol2.so: pyMol2.cpp pyMol2.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMol2.cpp
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMol2.h
	@echo 'Building Mol2 dinamic library.'
	$(CC) $< $(CC_FLAGS) pyMol2.cpp -o pyMol2.so $(LD_FLAGS)
	@echo ' '
	
install:
	
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMol2.cpp
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyMol2.h
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyWRITER.cpp
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyWRITER.h
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyGrid.h
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyGrid.cpp
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/iMcLiBELa.h
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyCOORD_MC.cpp
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyCOORD_MC.h
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyRAND.h
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyFindHB.cpp
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyFindHB.h
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyEnergy2.cpp
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/src/pyEnergy2.h
	
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/rec.mol2.gz
	@ wget https://raw.githubusercontent.com/alessandronascimento/pyLiBELa/main/test/lig.mol2.gz
	

	@ mv *.h *.cpp src/
	@ mv *gz test/

uninstall:
	@ rm -rf src
	@ rm -rf test
	
objFolder:
	@ mkdir -p obj
	@ mkdir src
	@ mkdir test
	
	
clean:
	@ $(RM) ./obj/*.so $(PROJ_NAME) *~
	@ rmdir obj

.PHONY: all clean
