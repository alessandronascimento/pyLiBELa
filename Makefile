#Name of the project
PROJ_NAME=pyLiBELa

#Every .cpp file
CPP_SRC=$(wildcard ./src/*.cpp)

#Every .h file
H_SRC=$(wildcard ./src/*.h)

#Dinamic libraries
OBJ=$(subst .cpp,.so,$(subst src,obj,$(CPP_SRC)))


#Compiler
CC=g++

#Flags for compiler
CC_FLAGS=-fPIC                       \
	 -shared                     \
	 -I/usr/include/python3.10   \
	 -I/usr/include/openbabel3   \
	 -I/usr/include/eigen3       \
	 -DBUILD=0                   \
	 -DHAVE_EIGEN                
	    
#Linkage
LD_FLAGS=-L/usr/lib/x86_64-linux-gnu \
	 -lboost_python310           \
	 -lz                         \
	 -lopenbabel                  


#Command to clean target
RM = rm -rf

#
#Compilation and linking
#
all: objFolder $(PROJ_NAME)

$(PROJ_NAME): $(OBJ)

./obj/%.so: ./src/%.cpp ./src/%.h
	@echo 'Building dinamic libraries using G++ compiler: $<'
	$(CC) $< $(CC_FLAGS) -o $@ $(LD_FLAGS)
	@echo ' '

objFolder:
	@ mkdir -p obj
	
	
clean:
	@ $(RM) ./obj/*.so $(PROJ_NAME) *~
	@ rmdir obj

.PHONY: all clean