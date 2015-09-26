#pragma once

#ifdef MAP_MODE
#include <sys/mman.h> // for mmap() / munmap()
#include <assert.h>
#define MEMORY_ALIGNMENT  4096
#define ALIGN_UP(x,size) ( ((size_t)x+(size-1))&(~(size-1)) )

inline void
AllocateHostMemory(bool bPinGenericMemory, void **pp_a, void **ppAligned_a, int nbytes)
{
	if (bPinGenericMemory)
	{
		// allocate a generic page-aligned chunk of system memory
		//printf("> mmap() allocating %4.2f Mbytes (generic page-aligned system memory)\n", (float)nbytes/1048576.0f);
		*pp_a = (void *)mmap(NULL, (nbytes + MEMORY_ALIGNMENT), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);

		*ppAligned_a = (void *)ALIGN_UP(*pp_a, MEMORY_ALIGNMENT);

		//printf("> cudaHostRegister() registering %4.2f Mbytes of generic allocated system memory\n", (float)nbytes/1048576.0f);
		// pin allocate memory
		cudaHostRegister(*ppAligned_a, nbytes, cudaHostRegisterMapped);
	}
	else
	{
		//printf("> cudaMallocHost() allocating %4.2f Mbytes of system memory\n", (float)nbytes/1048576.0f);
		// allocate host memory (pinned is required for achieve asynchronicity)
		cudaMallocHost((void **)pp_a, nbytes);
		*ppAligned_a = *pp_a;
	}
}

inline void
FreeHostMemory(bool bPinGenericMemory, void **pp_a, void **ppAligned_a, int nbytes)
{
	// CUDA 4.0 support pinning of generic host memory
	if (bPinGenericMemory)
	{
		// unpin and delete host memory
		cudaHostUnregister(*ppAligned_a);

		munmap(*pp_a, nbytes);
	}
	else
	{
		cudaFreeHost(*pp_a);
	}
}
#endif