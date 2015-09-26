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

#define MINLOG 1e-6 
#define MAXLOG 512

#define CODE1
//#define CODE2

#define LINUX

//#define TEST_PERF



#define ThreadpBlock 256
#define ITERATE_TIME 16	//Iteration times
#define STREAM_NUM 3	//Cuda_Streams

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

typedef long int __int64;

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

typedef float buf_info_col[BLOCK_SIZE];
typedef struct
{
	buf_info_col info[MAX_DEG_COL];
}BUF_INFO_COL;
typedef float buf_info_ch[BLOCK_SIZE];

typedef float info_ch[COL_LENGTH];


//Device Matrix 
__constant__ matrix_check_node d_matrix_node_c[BLOCK_NUM_ROW];
__constant__ matrix_variable_node d_matrix_node_v[BLOCK_NUM_COL];
