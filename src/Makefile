decode: kernel.o CPU_decode.o
	/usr/local/cuda-7.0/bin/nvcc -ccbin g++ -o decode kernel.o CPU_decode.o

kernel.o: kernel.cu totalDefine.h linux.h randn.h cuda_helper.cuh GPU_decode_10240.cuh GPU_decode_7168.cuh
	/usr/local/cuda-7.0/bin/nvcc -ccbin g++ -o kernel.o -arch=compute_20 -code=compute_20,sm_21 --machine 64 --use_fast_math --maxrregcount=32 -c kernel.cu

CPU_decode.o: CPU_decode.cpp
	/usr/local/cuda-7.0/bin/nvcc -ccbin g++ -c CPU_decode.cpp

clean:
	rm decode kernel.o CPU_decode.o
