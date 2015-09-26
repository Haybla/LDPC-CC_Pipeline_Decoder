#include "totalDefine.h"
#include "CPU_decode.h"
#include "randn.h"
#include "cuda_randn.cuh"

#ifdef CODE1
#include "GPU_decode_10240.cuh"
#endif
#ifdef CODE2
#include "GPU_decode_7168.cuh"
#endif

#ifdef MAP_MODE
#include "linux.h"
#endif

//Host Matrix
matrix_check_node h_matrix_node_c[BLOCK_NUM_ROW];
matrix_variable_node h_matrix_node_v[BLOCK_NUM_COL];


int main()
{
	/*****************************************************************/
	/********************GPU Device Initialization********************/
	/*****************************************************************/
	cudaDeviceReset();

#ifdef MAP_MODE
	cudaSetDeviceFlags(cudaDeviceMapHost);
#endif

	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);


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

#ifdef TESTDB_MODE
	__int64 TIME_COUNT = 1e9/STREAM_NUM/STREAM_COUNT;
	__int64 RELAY = 1e5;
	int MAX_FE = 500;
#else
	__int64 TIME_COUNT = 1e3;
#endif

#ifdef TESTDB_MODE
	float DB_start = 1.6;
	float DB_end = 1.7;
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
#ifdef BUFFER_MODE
		BUF_INFO_COL *buf_d_info_col_2_row[STREAM_NUM * SUB_NUM];
		BUF_INFO_COL *buf_h_info_col_2_row[STREAM_NUM * SUB_NUM];
		buf_info_ch *buf_d_channel_info[STREAM_NUM];
		buf_info_ch *buf_h_channel_info[STREAM_NUM];
#endif
#if defined(BUFFER_MODE) && defined(MAP_MODE)
		BUF_INFO_COL *buf_h_col[STREAM_NUM * SUB_NUM];
		buf_info_ch *buf_h_ch[STREAM_NUM * SUB_NUM];
		int *th_decoded_w[STREAM_NUM];
#endif

		//Malloc Host Memory
		for (int i = 0; i < (ITERATE_TIME * SUB_NUM); i++)
		{
			cudaHostAlloc(&h_decoded_word[i], sizeof(int)*BLOCK_SIZE, cudaHostAllocDefault);
		}
		cudaHostAlloc(&h_channel_info, sizeof(float)*COL_LENGTH, cudaHostAllocDefault);
		cudaHostAlloc(&h_info_col_2_row, sizeof(INFO_COL), cudaHostAllocDefault);
		cudaHostAlloc(&h_info_row_2_col, sizeof(INFO_ROW), cudaHostAllocDefault);

		//Malloc Device Memory and Host Temp Memory
		for (int i = 0; i < STREAM_NUM; i++)
		{
			cudaMalloc(&d_channel_info[i], sizeof(info_ch)*STREAM_COUNT);
			cudaMalloc(&d_decoded_word[i], sizeof(int)*BLOCK_SIZE*STREAM_COUNT);
			cudaMalloc(&d_info_col_2_row[i], sizeof(INFO_COL)*STREAM_COUNT);
			cudaMalloc(&d_info_row_2_col[i], sizeof(INFO_ROW)*STREAM_COUNT);

			cudaHostAlloc(&th_info_col_2_row[i], sizeof(INFO_COL)*STREAM_COUNT, cudaHostAllocDefault);
#ifndef MAP_MODE
			cudaHostAlloc(&th_decoded_word[i], sizeof(int)*BLOCK_SIZE*STREAM_COUNT, cudaHostAllocDefault);
#endif
		}

		//Malloc Buffers
#ifdef BUFFER_MODE
#ifdef MAP_MODE
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
			cudaHostAlloc(&buf_h_info_col_2_row[i * SUB_NUM + j], sizeof(BUF_INFO_COL)*STREAM_COUNT, cudaHostAllocDefault);
		}

		for (int i = 0; i < STREAM_NUM; i++)
		{
			cudaHostAlloc(&buf_h_channel_info[i], sizeof(buf_info_ch)*STREAM_COUNT * SUB_NUM, cudaHostAllocDefault);
		}
#endif
		for (int i = 0; i < STREAM_NUM; i++)
		for (int j = 0; j < SUB_NUM; j++)
		{
			cudaMalloc(&buf_d_info_col_2_row[i * SUB_NUM + j], sizeof(BUF_INFO_COL)*STREAM_COUNT);
		}

		for (int i = 0; i < STREAM_NUM; i++)
		{
			cudaMalloc(&buf_d_channel_info[i], sizeof(buf_info_ch)*STREAM_COUNT * SUB_NUM);
		}
#endif

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

#ifdef TESTDB_MODE
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
				cudaMemcpy(&d_channel_info[i][j], h_channel_info, sizeof(info_ch), cudaMemcpyHostToDevice);
				cudaMemcpy(&d_info_col_2_row[i][j], h_info_col_2_row, sizeof(INFO_COL), cudaMemcpyHostToDevice); 
				cudaMemcpy(&d_info_row_2_col[i][j], h_info_row_2_col, sizeof(INFO_ROW), cudaMemcpyHostToDevice); 
			}
		}

		//Matrix Copy Host to Device, Constant Memory
		cudaMemcpyToSymbol(d_matrix_node_c, h_matrix_node_c, sizeof(matrix_check_node)* BLOCK_NUM_ROW); 
		cudaMemcpyToSymbol(d_matrix_node_v, h_matrix_node_v, sizeof(matrix_variable_node)* BLOCK_NUM_COL);

		//Buffers Initialization
#ifdef BUFFER_MODE
		for (int str_count = 0; str_count < STREAM_NUM; str_count++)
		for (int bl = 0; bl < SUB_NUM; bl++)
		{
			for (int i = 0; i < STREAM_COUNT; i++)
			{
				memcpy(&buf_h_channel_info[str_count][bl * STREAM_COUNT + i][0], &h_channel_info[(((time_count) % (ITERATE_TIME * SUB_NUM)) - SUB_NUM + bl) * BLOCK_SIZE], sizeof(float)*BLOCK_SIZE);
			}
			cudaMemcpy(&buf_d_channel_info[str_count][bl * STREAM_COUNT], &buf_h_channel_info[str_count][bl * STREAM_COUNT], sizeof(buf_info_ch)* STREAM_COUNT, cudaMemcpyHostToDevice);
		}
#endif

		/*****************************************************************/
		/**************************GPU Decoding***************************/
		/*****************************************************************/
		//Kernel Dimension Setting
		dim3 block(ThreadpBlock);
		dim3 grid(ITERATE_TIME, STREAM_COUNT);

		//Streams Creation
		cudaStream_t *str = (cudaStream_t *)malloc(STREAM_NUM * sizeof(cudaStream_t));
		for (int i = 0; i < STREAM_NUM; i++)
			cudaStreamCreate(&str[i]);

#ifdef TESTDB_MODE
		printf("time_count = %13d", time_count);
		bit_error = 0; block_error = 0; block_num = 0;		//only count BER/FER simulation on GPU
#endif

		//GPU Decoding Starting
		
		//GPU Timer
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);
		while (time_count < TIME_COUNT
#ifdef TESTDB_MODE
			&& block_error < MAX_FE
#endif
			)
		{
			int bl = (time_count + 1) % SUB_NUM;

			//Generate buffers
#if defined(BUFFER_MODE) && !defined(TESTDB_MODE)
				for (int str_count = 0; str_count < STREAM_NUM; str_count++)
				for (int i = 0; i < STREAM_COUNT; i++)
				{
					memcpy(&buf_h_channel_info[str_count][bl*STREAM_COUNT + i][0], &h_channel_info[((time_count) % (ITERATE_TIME * SUB_NUM)) * BLOCK_SIZE], sizeof(float)*BLOCK_SIZE);
				}
#endif
#if defined(BUFFER_MODE) && defined(TESTDB_MODE)
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
				cudaMemcpyAsync(&buf_d_channel_info[str_count][bl*STREAM_COUNT], &buf_h_channel_info[str_count][bl*STREAM_COUNT], sizeof(buf_info_ch)*STREAM_COUNT, cudaMemcpyHostToDevice, str[str_count]);
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
				cudaMemcpyAsync(th_decoded_word[str_count], d_decoded_word[str_count], sizeof(int)*BLOCK_SIZE*STREAM_COUNT, cudaMemcpyDeviceToHost, str[str_count]);
			}

			//Cuda Streams Synchronizing
			cudaDeviceSynchronize();


			//Compute number of errors in test mode
#ifdef TESTDB_MODE
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
					fprintf(f1, "SNR=%.3f, BER=%.3e, FER=%.3e, BLOCK_NUM=%d, BLOCK_ERR=%d\n", DB, BER, FER, block_num, block_error);
					fclose(f1);
					f1 = NULL;
				}
			}
#endif

#if !defined(TESTDB_MODE) && !defined(TESTSTREAM_MODE)
			time_count++;
#else
			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b");
			printf("%13d", time_count);
			time_count = time_count + 1;
#endif
		}
		//GPU Decoding Ending

		//GPU Timer
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&testTime, start, stop);
		totalTime += testTime;

		//Compute BER and FER in test mode
#ifdef TESTDB_MODE
		float BER = (float)(bit_error) / (block_num * BLOCK_SIZE);
		float FER = (float)block_error / block_num;
		printf("\nBER is %.3e \nFER is %.3e \n", BER, FER);
#endif

		//Compute decoding time and Throughput
		printf("Time is %.3f ms\n", totalTime);
		printf("Thoughput is %.3f Mbps\n", (float)BLOCK_SIZE*(time_count-ITERATE_TIME*SUB_NUM)*STREAM_NUM*STREAM_COUNT / totalTime / 1000);
		printf("-------------------------------\n");

#ifdef TESTDB_MODE
		f1 = fopen("BER&FER.txt", "a+");
		fprintf(f1, "SNR=%.3f, BER=%.3e, FER=%.3e, BLOCK_NUM=%d, BLOCK_ERR=%d\n", DB, BER, FER, block_num, block_error);
		fclose(f1);
		f1 = NULL;
#endif

		/*****************************************************************/
		/*************************Memory Releasing************************/
		/*****************************************************************/
		//Free Host Memory
		for (int i = 0; i < (ITERATE_TIME * SUB_NUM); i++)
		{
			cudaFreeHost(h_decoded_word[i]);
		}
		cudaFreeHost(h_channel_info);
		cudaFreeHost(h_info_col_2_row);
		cudaFreeHost(h_info_row_2_col);

		//Free Device Memory and Temp Memory
		for (int i = 0; i < STREAM_NUM; i++)
		{
			cudaFree(d_channel_info[i]);
			cudaFree(d_decoded_word[i]);
			cudaFree(d_info_col_2_row[i]);
			cudaFree(d_info_row_2_col[i]);

			cudaFreeHost(th_info_col_2_row[i]);
#ifndef MAP_MODE
			cudaFreeHost(th_decoded_word[i]);
#endif
		}

		//Free Buffers
#ifdef BUFFER_MODE
#ifdef MAP_MODE
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
				cudaFreeHost(buf_h_info_col_2_row[i * SUB_NUM + j]);
			cudaFreeHost(buf_h_channel_info[i]);
		}
#endif
		for (int i = 0; i < STREAM_NUM; i++)
		{
			for (int j = 0; j < SUB_NUM; j++)
				cudaFree(buf_d_info_col_2_row[i * SUB_NUM + j]);
			cudaFree(buf_d_channel_info[i]);
		}
#endif

		//Cuda Stream Destroy
		for (int i = 0; i < STREAM_NUM; i++)
			cudaStreamDestroy(str[i]);

		//Cuda Event Destroy
		cudaEventDestroy(start);
		cudaEventDestroy(stop);

#ifdef TESTDB_MODE
	}
#endif

		//GPU Device Reset
		cudaDeviceReset();

		//Exit
		printf("Test passed\n");
		exit(EXIT_SUCCESS);

		return 0;
}