CUDA_FILES= *.cu

CPP_FILES= *.cpp

OBJS= *.o

EXEC= exec

CUDA_FLAGS= -L/usr/local/cuda/lib64 -lcudart

all: cuda cpp link

cuda: 
	@ nvcc -c $(CUDA_FILES)

cpp: 
	@ g++ -c $(CPP_FILES)

link: 
	@ nvcc -o $(EXEC) $(OBJS) -lboost_python310 

boost:
	@ g++ -Wall -Wl,--no-undefined -Wextra -fPIC -shared -I/usr/include/python3.10 $(CPP_FILES) -o test.so -lboost_python310


test: 
	@ g++ -c -I /usr/include/python3.10 testSource.cpp

run: 
	@./$(EXEC)

clear: 
	@ rm -f *.o
	@ rm -f $(EXEC)

