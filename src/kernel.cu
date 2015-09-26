/*
Copyright (c) 2014-2015 Mokky and Haybla. All rights reserved.

This file is part of LDPC-CC_Pipeline_Decoder. Original Codes can 
be found at <https://github.com/Haybla>.

LDPC-CC_Pipeline_Decoder is free software: you can redistribute it 
and/or modify it under the terms of the GNU General Public License 
as published by the Free Software Foundation, either version 3 of 
the License, or any later version.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "totalDefine.h"
#include "CPU_decode.h"
#include "randn.h"
#include "cuda_helper.cuh"

#ifdef CODE1
#include "GPU_decode_10240.cuh"
#endif
#ifdef CODE2
#include "GPU_decode_7168.cuh"
#endif

#ifdef LINUX
#include "linux.h"
#endif

#ifdef CODE1
#define PRINTPARM printf("CCSDS(%d,%d): ",4096, 10240);
#else
#define PRINTPARM printf("CCSDS(%d,%d): ",4096, 7168);
#endif

//Host Matrix
matrix_check_node h_matrix_node_c[BLOCK_NUM_ROW];
matrix_variable_node h_matrix_node_v[BLOCK_NUM_COL];


int main()
{
	/*****************************************************************/
	/********************GPU Device Initialization********************/
	/*****************************************************************/
	checkCudaErrors(cudaDeviceReset());

#ifdef LINUX
	checkCudaErrors(cudaSetDeviceFlags(cudaDeviceMapHost));
#endif

	checkCudaErrors(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));


	/*****************************************************************/
	/********************Simulation Initialization********************/
	/*****************************************************************/
	//Matrix Initialization
	fun_matrix();

	//Random Seed Initialization
	srand((unsigned)time(NULL));

	//GPU Timer Initialization
	cudaEvent_t start, stop;

	float totalTime = 0;
	float testTime = 0;

	//Simulation Parameters Setting: SNR, Code_Streams
	float DB = 1.2;
	int STREAM_COUNT = 15;

#ifdef TEST_PERF
	__int64 TIME_COUNT = 1e9/STREAM_NUM/STREAM_COUNT;
	__int64 RELAY = 1e5;
	int MAX_FE = 500;
#else
	__int64 TIME_COUNT = 1e4;
#endif

#ifdef TEST_PERF
	float DB_start = 1.1;
	float DB_end = 1.2;
	float DB_step = 0.05;

	static FILE *f1 = NULL;
	f1 = fopen("BER&FER.txt", "a+");
#ifdef CODE1
	fprintf(f1, "CODE=%d, DB_START=%1.3f, DB_END=%1.3f, DB_STEP=%1.3f, ITR=%d\n", 1, DB_start, DB_end, DB_step, ITERATE_TIME);
#endif
#ifdef CODE2
	fprintf(f1, "CODE=%d, DB_START=%1.3f, DB_END=%1.3f, DB_STEP=%1.3f, ITR=%d\n", 2, DB_start, DB_end, DB_step, ITERATE_TIME);
#endif
	fclose(f1);
	f1 = NULL;

	for (DB = DB_start; DB <= DB_end; DB += DB_step)
	{
#endif

		/*****************************************************************/
		/**********************Memory Initialization**********************/
		/*****************************************************************/
		//Host Memory
		float *h_channel_info; 
		int *h_decoded_word[ITERATE_TIME * SUB_NUM]; 
		INFO_COL *h_info_col_2_row;
		INFO_ROW *h_info_row_2_col;

		//Device Memory
		info_ch *d_channel_info[STREAM_NUM]; 
		INFO_COL *d_info_col_2_row[STREAM_NUM];
		INFO_ROW *d_info_row_2_col[STREAM_NUM];
		int *d_decoded_word[STREAM_NUM];

		//Host Temp Memory
		INFO_COL *th_info_col_2_row[STREAM_NUM];
		int *th_decoded_word[STREAM_NUM];

		//Buffers
		BUF_INFO_COL *buf_d_info_col_2_row[STREAM_NUM * SUB_NUM];
		BUF_INFO_COL *buf_h_info_col_2_row[STREAM_NUM * SUB_NUM];
		buf_info_ch *buf_d_channel_info[STREAM_NUM];
		buf_info_ch *buf_h_channel_info[STREAM_NUM];

#ifdef LINUX
		BUF_INFO_COL *buf_h_col[STREAM_NUM * SUB_NUM];
		buf_info_ch *buf_h_ch[STREAM_NUM * SUB_NUM];
		int *th_decoded_w[STREAM_NUM];
#endif

		//Malloc Host Memory
		for (int i = 0; i < (ITERATE_TIME * SUB_NUM); i++)
		{
			checkCudaErrors(cudaHostAlloc(&h_decoded_word[i], sizeof(int)*BLOCK_SIZE, cudaHostAllocDefault));
		}
		checkCudaErrors(cudaHostAlloc(&h_channel_info, sizeof(float)*COL_LENGTH, cudaHostAllocDefault));
		checkCudaErrors(cudaHostAlloc(&h_info_col_2_row, sizeof(INFO_COL), cudaHostAllocDefault));
		checkCudaErrors(cudaHostAlloc(&h_info_row_2_col, sizeof(INFO_ROW), cudaHostAllocDefault));

		//Malloc Device Memory and Host Temp Memory
		for (int i = 0; i < STREAM_NUM; i++)
		{
			checkCudaErrors(cudaMalloc(&d_channel_info[i], sizeof(info_ch)*STREAM_COUNT));
			checkCudaErrors(cudaMalloc(&d_decoded_word[i], sizeof(int)*BLOCK_SIZE*STREAM_COUNT));
			checkCudaErrors(cudaMalloc(&d_info_col_2_row[i], sizeof(INFO_COL)*STREAM_COUNT));
			checkCudaErrors(cudaMalloc(&d_info_row_2_col[i], sizeof(INFO_ROW)*STREAM_COUNT));

			checkCudaErrors(cudaHostAlloc(&th_info_col_2_row[i], sizeof(INFO_COL)*STREAM_COUNT, cudaHostAllocDefault));
#ifndef LINUX
			checkCudaErrors(cudaHostAlloc(&th_decoded_word[i], sizeof(int)*BLOCK_SIZE*STREAM_COUNT, cudaHostAllocDefault));
#endif
		}

		//Malloc Buffers
#ifdef LINUX
		for (int i = 0; i < STREAM_NUM; i++)
		for (int j = 0; j < SUB_NUM; j++)
		{
			AllocateHostMemory(1, (void **)&buf_h_col[i * SUB_NUM + j], (void **)&buf_h_info_col_2_row[i * SUB_NUM + j], sizeof(BUF_INFO_COL)*STREAM_COUNT);
		}

		for (int i = 0; i < STREAM_NUM; i++)
		{
			AllocateHostMemory(1, (void **)&buf_h_ch[i], (void **)&buf_h_channel_info[i], sizeof(buf_info_ch)*STREAM_COUNT * 4);
		}

		for (int i = 0; i < STREAM_NUM; i++)
		{
			AllocateHostMemory(1, (void **)&th_decoded_w[i], (void **)&th_decoded_word[i], sizeof(int)*BLOCK_SIZE*STREAM_COUNT);
		}
#else
		for (int i = 0; i < STREAM_NUM; i++)
		for (int j = 0; j < SUB_NUM; j++)
		{
			checkCudaErrors(cudaHostAlloc(&buf_h_info_col_2_row[i * SUB_NUM + j], sizeof(BUF_INFO_COL)*STREAM_COUNT, cudaHostAllocDefault));
		}

		for (int i = 0; i < STREAM_NUM; i++)
		{
			checkCudaErrors(cudaHostAlloc(&buf_h_channel_info[i], sizeof(buf_info_ch)*STREAM_COUNT * SUB_NUM, cudaHostAllocDefault));
		}
#endif
		for (int i = 0; i < STREAM_NUM; i++)
		for (int j = 0; j < SUB_NUM; j++)
		{
			checkCudaErrors(cudaMalloc(&buf_d_info_col_2_row[i * SUB_NUM + j], sizeof(BUF_INFO_COL)*STREAM_COUNT));
		}

		for (int i = 0; i < STREAM_NUM; i++)
		{
			checkCudaErrors(cudaMalloc(&buf_d_channel_info[i], sizeof(buf_info_ch)*STREAM_COUNT * SUB_NUM));
		}


		/*****************************************************************/
		/***********************Simulation Starting***********************/
		/*****************************************************************/
		//Change SNR to Sigma^2
		float no = NoCal(DB);

		//Error Counters Initialization
		int err = 0;
		__int64 block_error = 0;
		__int64 block_num = 0;
		__int64 bit_error = 0;
		__int64 time_count = 0;
		totalTime = 0;

#ifdef TEST_PERF
		printf("-----dB is %1.2f Testing Now-----\n", DB);
#endif


		/*****************************************************************/
		/**************************CPU Decoding***************************/
		/*****************************************************************/
		//The first ITERATE_TIME*SUB_NUM-1 is decoding on CPU
		while (time_count < (ITERATE_TIME * SUB_NUM - 1))
		{
			V_rand(no, &h_channel_info[((time_count) % (ITERATE_TIME * SUB_NUM)) * BLOCK_SIZE]);

			for (int i = 0; i < MAX_DEG_COL; i++)
			{
				memcpy(&h_info_col_2_row->info[i][((time_count) % (ITERATE_TIME * SUB_NUM)) * BLOCK_SIZE], &h_channel_info[((time_count) % (ITERATE_TIME * SUB_NUM)) * BLOCK_SIZE], sizeof(float)*BLOCK_SIZE);
			}

			update(time_count, h_channel_info, h_decoded_word, h_info_col_2_row, h_info_row_2_col);

			if (time_count >= (ITERATE_TIME * SUB_NUM - 1))
			{
				err = countError(h_decoded_word[(time_count - (ITERATE_TIME * SUB_NUM - 1)) % (ITERATE_TIME * SUB_NUM)]);
				bit_error += err;
				if (err > 0)
				{
					block_error++;
				}
				block_num++;
			}

			time_count++;
		}
		//CPU Decoding Ending


		//Memory Copy Host to Device, Global Memory
		for (int i = 0; i < STREAM_NUM; i++)
		{
			for (int j = 0; j < STREAM_COUNT; j++)
			{
				checkCudaErrors(cudaMemcpy(&d_channel_info[i][j], h_channel_info, sizeof(info_ch), cudaMemcpyHostToDevice));
				checkCudaErrors(cudaMemcpy(&d_info_col_2_row[i][j], h_info_col_2_row, sizeof(INFO_COL), cudaMemcpyHostToDevice));
				checkCudaErrors(cudaMemcpy(&d_info_row_2_col[i][j], h_info_row_2_col, sizeof(INFO_ROW), cudaMemcpyHostToDevice));
			}
		}

		//Matrix Copy Host to Device, Constant Memory
		checkCudaErrors(cudaMemcpyToSymbol(d_matrix_node_c, h_matrix_node_c, sizeof(matrix_check_node)* BLOCK_NUM_ROW));
		checkCudaErrors(cudaMemcpyToSymbol(d_matrix_node_v, h_matrix_node_v, sizeof(matrix_variable_node)* BLOCK_NUM_COL));

		//Buffers Initialization
		for (int str_count = 0; str_count < STREAM_NUM; str_count++)
		for (int bl = 0; bl < SUB_NUM; bl++)
		{
			for (int i = 0; i < STREAM_COUNT; i++)
			{
				memcpy(&buf_h_channel_info[str_count][bl * STREAM_COUNT + i][0], &h_channel_info[(((time_count) % (ITERATE_TIME * SUB_NUM)) - SUB_NUM + bl) * BLOCK_SIZE], sizeof(float)*BLOCK_SIZE);
			}
			checkCudaErrors(cudaMemcpy(&buf_d_channel_info[str_count][bl * STREAM_COUNT], &buf_h_channel_info[str_count][bl * STREAM_COUNT], sizeof(buf_info_ch)* STREAM_COUNT, cudaMemcpyHostToDevice));
		}


		/*****************************************************************/
		/**************************GPU Decoding***************************/
		/*****************************************************************/
		//Kernel Dimension Setting
		dim3 block(ThreadpBlock);
		dim3 grid(ITERATE_TIME, STREAM_COUNT);

		//Streams Creation
		cudaStream_t *str = (cudaStream_t *)malloc(STREAM_NUM * sizeof(cudaStream_t));
		for (int i = 0; i < STREAM_NUM; i++)
			checkCudaErrors(cudaStreamCreate(&str[i]));

#ifdef TEST_PERF
		printf("time_count = %13ld", time_count);
		bit_error = 0; block_error = 0; block_num = 0;		//only count BER/FER simulation on GPU
#endif

		//GPU Decoding Starting
		
		//GPU Timer
		checkCudaErrors(cudaEventCreate(&start));
		checkCudaErrors(cudaEventCreate(&stop));

		checkCudaErrors(cudaEventRecord(start, 0));
		while (time_count < TIME_COUNT
#ifdef TEST_PERF
			&& block_error < MAX_FE
#endif
			)
		{
			int bl = (time_count + 1) % SUB_NUM;

			//Generate buffers
#ifndef TEST_PERF
				for (int str_count = 0; str_count < STREAM_NUM; str_count++)
				for (int i = 0; i < STREAM_COUNT; i++)
				{
					memcpy(&buf_h_channel_info[str_count][bl*STREAM_COUNT + i][0], &h_channel_info[((time_count) % (ITERATE_TIME * SUB_NUM)) * BLOCK_SIZE], sizeof(float)*BLOCK_SIZE);
				}
#else
				for (int str_count = 0; str_count < STREAM_NUM; str_count++)
				for (int i = 0; i < STREAM_COUNT; i++)
				{
					//Generate all-zero codeword with Guass-Noise, length of BLOCK_SIZE
					V_rand(no, &h_channel_info[((time_count) % (ITERATE_TIME * SUB_NUM)) * BLOCK_SIZE]);
					memcpy(&buf_h_channel_info[str_count][bl*STREAM_COUNT + i][0], &h_channel_info[((time_count) % (ITERATE_TIME * SUB_NUM)) * BLOCK_SIZE], sizeof(float)*BLOCK_SIZE);
				}
#endif


			//H2D
			for (int str_count = 0; str_count < STREAM_NUM; str_count++)
			{
				checkCudaErrors(cudaMemcpyAsync(&buf_d_channel_info[str_count][bl*STREAM_COUNT], &buf_h_channel_info[str_count][bl*STREAM_COUNT], sizeof(buf_info_ch)*STREAM_COUNT, cudaMemcpyHostToDevice, str[str_count]));
			}

			//Kernel Execution
			for (int str_count = 0; str_count < STREAM_NUM; str_count++)
			{
				switch (bl)
				{
				case 0:
				{
						  dUpdate1 << <grid, block, 0, str[str_count] >> >(time_count, d_channel_info[str_count], d_decoded_word[str_count], d_info_col_2_row[str_count], d_info_row_2_col[str_count]
							  , buf_d_channel_info[str_count], buf_d_info_col_2_row[str_count * 4], STREAM_COUNT
							  );
						  break;
				}
				case 1:
				{
						  dUpdate2 << <grid, block, 0, str[str_count] >> >(time_count, d_channel_info[str_count], d_decoded_word[str_count], d_info_col_2_row[str_count], d_info_row_2_col[str_count]
							  , buf_d_channel_info[str_count], buf_d_info_col_2_row[str_count * 4], STREAM_COUNT
							  );
						  break;
				}
				case 2:
				{
						  dUpdate3 << <grid, block, 0, str[str_count] >> >(time_count, d_channel_info[str_count], d_decoded_word[str_count], d_info_col_2_row[str_count], d_info_row_2_col[str_count]
							  , buf_d_channel_info[str_count], buf_d_info_col_2_row[str_count * 4], STREAM_COUNT
							  );
						  break;
				}
				case 3:
				{
						  dUpdate4 << <grid, block, 0, str[str_count] >> >(time_count, d_channel_info[str_count], d_decoded_word[str_count], d_info_col_2_row[str_count], d_info_row_2_col[str_count]
							  , buf_d_channel_info[str_count], buf_d_info_col_2_row[str_count * 4], STREAM_COUNT
							  );
						  break;
				}
				default:
					break;
				}
			}

			//D2H
			for (int str_count = 0; str_count < STREAM_NUM; str_count++)
			{
				checkCudaErrors(cudaMemcpyAsync(th_decoded_word[str_count], d_decoded_word[str_count], sizeof(int)*BLOCK_SIZE*STREAM_COUNT, cudaMemcpyDeviceToHost, str[str_count]));
			}

			//Cuda Streams Synchronizing
			checkCudaErrors(cudaDeviceSynchronize());


			//Compute number of errors in test mode
#ifdef TEST_PERF
			if (time_count > 1000)
			{
				for (int i = 0; i < STREAM_NUM; i++)
				for (int j = 0; j < STREAM_COUNT; j++)
				{
					err = countError(&th_decoded_word[i][BLOCK_SIZE*j]);
					bit_error += err;
					if (err > 0)
					{
						block_error++;
					}
				}

				block_num = block_num + STREAM_NUM*STREAM_COUNT;
				
				if (block_num % RELAY == 0)
				{
					float BER = (float)(bit_error) / (block_num * BLOCK_SIZE);
					float FER = (float)block_error / block_num;
					f1 = fopen("BER&FER.txt", "a+");
					fprintf(f1, "SNR=%.3f, BER=%.3e, FER=%.3e, BLOCK_NUM=%ld, BLOCK_ERR=%ld\n", DB, BER, FER, block_num, block_error);
					fclose(f1);
					f1 = NULL;
				}
			}
#endif

#ifndef TEST_PERF
			time_count++;
#else
			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b");
			printf("%13ld", time_count);
			time_count = time_count + 1;
#endif
		}
		//GPU Decoding Ending

		//GPU Timer
		checkCudaErrors(cudaEventRecord(stop, 0));
		checkCudaErrors(cudaEventSynchronize(stop));
		checkCudaErrors(cudaEventElapsedTime(&testTime, start, stop));
		totalTime += testTime;

		//Compute BER and FER in test mode
#ifdef TEST_PERF
		float BER = (float)(bit_error) / (block_num * BLOCK_SIZE);
		float FER = (float)block_error / block_num;
		printf("\n");
		PRINTPARM
		printf("BER is %.3e \n", BER);
		PRINTPARM
		printf("FER is %.3e \n", FER);
#endif

#ifndef TEST_PERF
		//Compute decoding time and Throughput
		PRINTPARM
		printf("Total Number of bits is %.3f Mb\n", (float)BLOCK_SIZE*(time_count-ITERATE_TIME*SUB_NUM)*STREAM_NUM*STREAM_COUNT/1024/1024);
		PRINTPARM		
		printf("Time is %.3f ms\n", totalTime);
		PRINTPARM
		printf("Thoughput is %.3f Mbps\n", (float)BLOCK_SIZE*(time_count-ITERATE_TIME*SUB_NUM)*STREAM_NUM*STREAM_COUNT / totalTime / 1000);
#endif		
		printf("-------------------------------\n");

#ifdef TEST_PERF
		f1 = fopen("BER&FER.txt", "a+");
		fprintf(f1, "SNR=%.3f, BER=%.3e, FER=%.3e, BLOCK_NUM=%ld, BLOCK_ERR=%ld\n", DB, BER, FER, block_num, block_error);
		fclose(f1);
		f1 = NULL;
#endif

		/*****************************************************************/
		/*************************Memory Releasing************************/
		/*****************************************************************/
		//Free Host Memory
		for (int i = 0; i < (ITERATE_TIME * SUB_NUM); i++)
		{
			checkCudaErrors(cudaFreeHost(h_decoded_word[i]));
		}
		checkCudaErrors(cudaFreeHost(h_channel_info));
		checkCudaErrors(cudaFreeHost(h_info_col_2_row));
		checkCudaErrors(cudaFreeHost(h_info_row_2_col));

		//Free Device Memory and Temp Memory
		for (int i = 0; i < STREAM_NUM; i++)
		{
			checkCudaErrors(cudaFree(d_channel_info[i]));
			checkCudaErrors(cudaFree(d_decoded_word[i]));
			checkCudaErrors(cudaFree(d_info_col_2_row[i]));
			checkCudaErrors(cudaFree(d_info_row_2_col[i]));

			checkCudaErrors(cudaFreeHost(th_info_col_2_row[i]));
#ifndef LINUX
			checkCudaErrors(cudaFreeHost(th_decoded_word[i]));
#endif
		}

		//Free Buffers
#ifdef LINUX
		for (int i = 0; i < STREAM_NUM; i++)
		{
			for (int j = 0; j < SUB_NUM; j++)
				FreeHostMemory(1, (void **)&buf_h_col[i * SUB_NUM + j], (void **)&buf_h_info_col_2_row[i * SUB_NUM + j], sizeof(BUF_INFO_COL)*STREAM_COUNT);
			FreeHostMemory(1, (void **)&buf_h_ch[i], (void **)&buf_h_channel_info[i], sizeof(buf_info_ch)*STREAM_COUNT * 4);
			FreeHostMemory(1, (void **)&th_decoded_w[i], (void **)&th_decoded_word[i], sizeof(int)*BLOCK_SIZE*STREAM_COUNT);
		}
#else
		for (int i = 0; i < STREAM_NUM; i++)
		{
			for (int j = 0; j < SUB_NUM; j++)
				checkCudaErrors(cudaFreeHost(buf_h_info_col_2_row[i * SUB_NUM + j]));
			checkCudaErrors(cudaFreeHost(buf_h_channel_info[i]));
		}
#endif
		for (int i = 0; i < STREAM_NUM; i++)
		{
			for (int j = 0; j < SUB_NUM; j++)
				checkCudaErrors(cudaFree(buf_d_info_col_2_row[i * SUB_NUM + j]));
			checkCudaErrors(cudaFree(buf_d_channel_info[i]));
		}

		//Cuda Stream Destroy
		for (int i = 0; i < STREAM_NUM; i++)
			checkCudaErrors(cudaStreamDestroy(str[i]));

		//Cuda Event Destroy
		checkCudaErrors(cudaEventDestroy(start));
		checkCudaErrors(cudaEventDestroy(stop));

#ifdef TEST_PERF
	}
#endif

		//GPU Device Reset
		checkCudaErrors(cudaDeviceReset());

		//Exit
		printf("Test passed\n");
		exit(EXIT_SUCCESS);

		//return 0;
}
