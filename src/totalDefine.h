#pragma once

#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "string.h"
#include "cuda_runtime.h"
#include "math.h"
#include "time.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
#include "cuda_device_runtime_api.h"
#include <curand.h>
#include <curand_kernel.h>

#define MINLOG 1e-6 
#define MAXLOG 512

//#define CODE1
#define CODE2

#define BUFFER_MODE
//#define CH_MODE
//#define MAP_MODE

#define COMPUTE_2x
//#define COMPUTE_3x

#define UNROLL

//#define TESTDB_MODE
//#define TESTSTREAM_MODE

#define ThreadpBlock 256
#define ITERATE_TIME 16	//Iteration times
#define STREAM_NUM 2	//Cuda_Streams

#ifdef CODE1
#define ROW_OFFSET 3
#define SUB_NUM 4
#define BLOCK_NUM_COL 20
#define BLOCK_NUM_ROW 12
#define ROW_LENGTH (ITERATE_TIME*6144)
#define COL_LENGTH (ITERATE_TIME*10240)
#define BLOCK_SIZE 2560
#define CHECK_SIZE 1536
#define MSG_SIZE 1024
#define QC_SIZE 512
#define MAX_DEG_COL 6
#define MAX_DEG_ROW 6
#endif

#ifdef CODE2
#define ROW_OFFSET 3
#define SUB_NUM 4
#define BLOCK_NUM_COL 28
#define BLOCK_NUM_ROW 12
#define ROW_LENGTH (ITERATE_TIME*3072)
#define COL_LENGTH (ITERATE_TIME*7168)
#define BLOCK_SIZE 1792
#define CHECK_SIZE 768
#define MSG_SIZE 1024
#define QC_SIZE 256
#define MAX_DEG_COL 6
#define MAX_DEG_ROW 10
#endif

typedef struct
{
	int pos;	//absolute position
	int cpos;	//relative position
}node_position;

typedef struct 
{
	int col; 
	int col_2_row; 
	int number;	 //offset
}check_node;

typedef struct
{
	int row_deg; 
	check_node order[MAX_DEG_ROW]; 
}matrix_check_node; 

typedef struct 
{
	int row; 
	int row_2_col; 
	int number;	 
}variable_node;

typedef struct
{
	int col_deg; 
	variable_node order[MAX_DEG_COL];
}matrix_variable_node;


typedef float info_col[COL_LENGTH];
typedef struct
{
	info_col info[6];
}INFO_COL; 

typedef float info_row[ROW_LENGTH];
typedef struct
{
	info_row info[MAX_DEG_ROW];
}INFO_ROW; 

#ifdef BUFFER_MODE
typedef float buf_info_col[BLOCK_SIZE];
typedef struct
{
	buf_info_col info[MAX_DEG_COL];
}BUF_INFO_COL;
typedef float buf_info_ch[BLOCK_SIZE];
#endif

typedef float info_ch[COL_LENGTH];


//Device Matrix 
__constant__ matrix_check_node d_matrix_node_c[BLOCK_NUM_ROW];
__constant__ matrix_variable_node d_matrix_node_v[BLOCK_NUM_COL];