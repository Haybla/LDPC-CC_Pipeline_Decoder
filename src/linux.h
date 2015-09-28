/*
Copyright (c) 2014-2015 Mokyy and Haybla. All rights reserved.

This file is part of LDPC-CC_Pipeline_Decoder. Original Codes can
be found at <https://github.com/Haybla>.

LDPC-CC_Pipeline_Decoder is free software: you can redistribute it
and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of
the License, or any later version.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#ifdef LINUX
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