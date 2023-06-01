CUDA_FILES= *.cu

CPP_FILES= *.cpp

OBJS= *.o

PYTHON_FILE = *.py

EXEC= cudaTest

CUDA_FLAGS= -L/usr/local/cuda/lib64 -lcudart
#NVCC= ~/miniconda3/bin/nvcc

BOOST_FLAGS = -lboost_python38 
# BOOST_FLAGS = -lboost_python39 -L ~/miniconda3/lib 
			  

#PYTHON_INCLUDE= -I ~/miniconda3/include/python3.9

PYTHON_INCLUDE= -I /usr/include/python3.8

all: cuda cpp linknvcc

python: 
	@python3.8 $(PYTHON_FILE)

cuda: 
	@ nvcc --compiler-options '-fPIC' -I /usr/local/cuda/include -c $(CUDA_FILES)

cpp: 
	@ g++ -fPIC -c $(PYTHON_INCLUDE) $(CPP_FILES)

linknvcc: 
	@ nvcc --shared $(OBJS) -o $(EXEC).so $(BOOST_FLAGS) $(CUDA_FLAGS)

linkg++:
	@ g++ -o $(EXEC) $(OBJS) $(BOOST_FLAGS) 

# boost:
# 	@ g++ -Wall -Wl,--no-undefined -Wextra -fPIC -shared -I/usr/include/python3.10 $(CPP_FILES) -o test.so -lboost_python310

run: 
	@./$(EXEC)

clear: 
	@ rm -f *.o
	@ rm -f *.so
	@ rm -f $(EXEC)

