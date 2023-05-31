CUDA_FILES= *.cu

CPP_FILES= *.cpp

OBJS= *.o

EXEC= exec

# CUDA_FLAGS= -L/usr/local/cuda/lib64 -lcudart

BOOST_FLAGS = -lboost_python38 

PYTHON_INCLUDE= -I /usr/include/python3.8

all: cuda cpp link

cuda: 
	@ nvcc -c $(CUDA_FILES)

cpp: 
	@ g++ -c $(PYTHON_INCLUDE) $(CPP_FILES)

linknvcc: 
	@ nvcc -o $(EXEC) $(OBJS) $(BOOST_FLAGS) 

linkg++:
	@ g++ -o $(EXEC) $(OBJS) $(BOOST_FLAGS) 

# boost:
# 	@ g++ -Wall -Wl,--no-undefined -Wextra -fPIC -shared -I/usr/include/python3.10 $(CPP_FILES) -o test.so -lboost_python310

run: 
	@./$(EXEC)

clear: 
	@ rm -f *.o
	@ rm -f $(EXEC)

