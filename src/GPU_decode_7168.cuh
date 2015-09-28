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

#include "totalDefine.h"

#ifdef CODE2

__noinline__ __device__ float dLtanh(float info)
{
	float y, x_abs, x_e;

	x_abs = (float)fabs(info);

	x_e = (float)__expf(-x_abs);
	y = (float)__logf((1 + x_e) / (1 - x_e));

	return (y);
};


__global__ void dUpdate1(__int64 time_count, info_ch *d_channel_info, int *d_decoded_word, INFO_COL *d_info_col_2_row, INFO_ROW *d_info_row_2_col
	, buf_info_ch *buf_d_channel_info, BUF_INFO_COL *buf_d_info_col_2_row, int STREAM_COUNT
	)
{
	int tid = threadIdx.x;
	int I = blockIdx.x;
	int Y = blockIdx.y;

	int row_offset = time_count * 768;
	int col_offset = (time_count - 3) * 1792;
	int row_id = 0;
	int col_id = 0;
	register int node_pos0_0, node_pos0_1, node_pos0_2, node_pos0_3, node_pos0_4, node_pos0_5, node_pos0_6, node_pos0_7, node_pos0_8, node_pos0_9;
	register int node_cpos0, node_cpos1, node_cpos2, node_cpos3, node_cpos4, node_cpos5, node_cpos6, node_cpos7, node_cpos8, node_cpos9;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	register float info0_0, info0_1, info0_2, info0_3, info0_4, info0_5, info0_6, info0_7, info0_8, info0_9; //info = channel col to row info
	register int info_symbol0, info_symbol1, info_symbol2, info_symbol3, info_symbol4, info_symbol5, info_symbol6, info_symbol7, info_symbol8, info_symbol9; //positive = 0, negative = 1
	register float info_sum0;
	register int info_sum_symbol0;
	register int symbol0;
	register float info_channel0;
	int node_pos_2 = 0;


	/***************************************************************/
	/***********************ROW_UPDATE******************************/
	/***************************************************************/

	/**************************ROW_0********************************/
	row_id = row_offset - I * 4 * 768 + 0 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 1 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 4 + (243 + tid) & 255];
		info0_2 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 5 + (42 + tid) & 255];
		info0_3 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 0 + (52 + tid) & 255];
		info0_4 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 2 + (0 + tid) & 255];
		info0_5 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 5 + (141 + tid) & 255];
		info0_6 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 1 + (70 + tid) & 255];
		info0_7 = buf_d_channel_info[Y][256 * 0 + (237 + tid) & 255];
		info0_8 = buf_d_channel_info[Y][256 * 1 + (77 + tid) & 255];
		info0_9 = buf_d_channel_info[Y][256 * 4 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 256 + ((d_matrix_node_c[block_id].order[3].number + tid) & 255)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 256 + ((d_matrix_node_c[block_id].order[4].number + tid) & 255)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 256 + ((d_matrix_node_c[block_id].order[5].number + tid) & 255)) % COL_LENGTH;
		node_pos0_6 = (offset + d_matrix_node_c[block_id].order[6].col * 256 + ((d_matrix_node_c[block_id].order[6].number + tid) & 255)) % COL_LENGTH;
		node_pos0_7 = (offset + d_matrix_node_c[block_id].order[7].col * 256 + ((d_matrix_node_c[block_id].order[7].number + tid) & 255)) % COL_LENGTH;
		node_pos0_8 = (offset + d_matrix_node_c[block_id].order[8].col * 256 + ((d_matrix_node_c[block_id].order[8].number + tid) & 255)) % COL_LENGTH;
		node_pos0_9 = (offset + d_matrix_node_c[block_id].order[9].col * 256 + ((d_matrix_node_c[block_id].order[9].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;
		node_cpos6 = (d_matrix_node_c[block_id].order[6].col_2_row) % COL_LENGTH;
		node_cpos7 = (d_matrix_node_c[block_id].order[7].col_2_row) % COL_LENGTH;
		node_cpos8 = (d_matrix_node_c[block_id].order[8].col_2_row) % COL_LENGTH;
		node_cpos9 = (d_matrix_node_c[block_id].order[9].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];
		info0_6 = d_info_col_2_row[Y].info[node_cpos6][node_pos0_6];
		info0_7 = d_info_col_2_row[Y].info[node_cpos7][node_pos0_7];
		info0_8 = d_info_col_2_row[Y].info[node_cpos8][node_pos0_8];
		info0_9 = d_info_col_2_row[Y].info[node_cpos9][node_pos0_9];
	}

	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info_symbol6 = ((unsigned int)(int)(info0_6 - 1)) >> 31;
	info0_6 = dLtanh(info0_6);
	info_sum_symbol0 ^= info_symbol6;
	info_sum0 += info0_6;

	info_symbol7 = ((unsigned int)(int)(info0_7 - 1)) >> 31;
	info0_7 = dLtanh(info0_7);
	info_sum_symbol0 ^= info_symbol7;
	info_sum0 += info0_7;

	info_symbol8 = ((unsigned int)(int)(info0_8 - 1)) >> 31;
	info0_8 = dLtanh(info0_8);
	info_sum_symbol0 ^= info_symbol8;
	info_sum0 += info0_8;

	info_symbol9 = ((unsigned int)(int)(info0_9 - 1)) >> 31;
	info0_9 = dLtanh(info0_9);
	info_sum_symbol0 ^= info_symbol9;
	info_sum0 += info0_9;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol5));

	info0_6 = info_sum0 - info0_6;
	info0_6 = dLtanh(info0_6);
	info0_6 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol6));

	info0_7 = info_sum0 - info0_7;
	info0_7 = dLtanh(info0_7);
	info0_7 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol7));

	info0_8 = info_sum0 - info0_8;
	info0_8 = dLtanh(info0_8);
	info0_8 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol8));

	info0_9 = info_sum0 - info0_9;
	info0_9 = dLtanh(info0_9);
	info0_9 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol9));



	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;
	d_info_row_2_col[Y].info[6][row_id%ROW_LENGTH] = info0_6;
	d_info_row_2_col[Y].info[7][row_id%ROW_LENGTH] = info0_7;
	d_info_row_2_col[Y].info[8][row_id%ROW_LENGTH] = info0_8;
	d_info_row_2_col[Y].info[9][row_id%ROW_LENGTH] = info0_9;



	/**************************ROW_1********************************/
	row_id = row_offset - I * 4 * 768 + 1 * 256 + tid;
	info_sum0 = 0; symbol0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 2 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 4 + (20 + tid) & 255];
		info0_2 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 5 + (179 + tid) & 255];
		info0_3 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 6 + (93 + tid) & 255];
		info0_4 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 3 + (0 + tid) & 255];
		info0_5 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 5 + (84 + tid) & 255];
		info0_6 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 6 + (206 + tid) & 255];
		info0_7 = buf_d_channel_info[Y][256 * 1 + (122 + tid) & 255];
		info0_8 = buf_d_channel_info[Y][256 * 2 + (67 + tid) & 255];
		info0_9 = buf_d_channel_info[Y][256 * 5 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 256 + ((d_matrix_node_c[block_id].order[3].number + tid) & 255)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 256 + ((d_matrix_node_c[block_id].order[4].number + tid) & 255)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 256 + ((d_matrix_node_c[block_id].order[5].number + tid) & 255)) % COL_LENGTH;
		node_pos0_6 = (offset + d_matrix_node_c[block_id].order[6].col * 256 + ((d_matrix_node_c[block_id].order[6].number + tid) & 255)) % COL_LENGTH;
		node_pos0_7 = (offset + d_matrix_node_c[block_id].order[7].col * 256 + ((d_matrix_node_c[block_id].order[7].number + tid) & 255)) % COL_LENGTH;
		node_pos0_8 = (offset + d_matrix_node_c[block_id].order[8].col * 256 + ((d_matrix_node_c[block_id].order[8].number + tid) & 255)) % COL_LENGTH;
		node_pos0_9 = (offset + d_matrix_node_c[block_id].order[9].col * 256 + ((d_matrix_node_c[block_id].order[9].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;
		node_cpos6 = (d_matrix_node_c[block_id].order[6].col_2_row) % COL_LENGTH;
		node_cpos7 = (d_matrix_node_c[block_id].order[7].col_2_row) % COL_LENGTH;
		node_cpos8 = (d_matrix_node_c[block_id].order[8].col_2_row) % COL_LENGTH;
		node_cpos9 = (d_matrix_node_c[block_id].order[9].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];
		info0_6 = d_info_col_2_row[Y].info[node_cpos6][node_pos0_6];
		info0_7 = d_info_col_2_row[Y].info[node_cpos7][node_pos0_7];
		info0_8 = d_info_col_2_row[Y].info[node_cpos8][node_pos0_8];
		info0_9 = d_info_col_2_row[Y].info[node_cpos9][node_pos0_9];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info_symbol6 = ((unsigned int)(int)(info0_6 - 1)) >> 31;
	info0_6 = dLtanh(info0_6);
	info_sum_symbol0 ^= info_symbol6;
	info_sum0 += info0_6;

	info_symbol7 = ((unsigned int)(int)(info0_7 - 1)) >> 31;
	info0_7 = dLtanh(info0_7);
	info_sum_symbol0 ^= info_symbol7;
	info_sum0 += info0_7;

	info_symbol8 = ((unsigned int)(int)(info0_8 - 1)) >> 31;
	info0_8 = dLtanh(info0_8);
	info_sum_symbol0 ^= info_symbol8;
	info_sum0 += info0_8;

	info_symbol9 = ((unsigned int)(int)(info0_9 - 1)) >> 31;
	info0_9 = dLtanh(info0_9);
	info_sum_symbol0 ^= info_symbol9;
	info_sum0 += info0_9;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol5));

	info0_6 = info_sum0 - info0_6;
	info0_6 = dLtanh(info0_6);
	info0_6 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol6));

	info0_7 = info_sum0 - info0_7;
	info0_7 = dLtanh(info0_7);
	info0_7 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol7));

	info0_8 = info_sum0 - info0_8;
	info0_8 = dLtanh(info0_8);
	info0_8 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol8));

	info0_9 = info_sum0 - info0_9;
	info0_9 = dLtanh(info0_9);
	info0_9 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol9));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;
	d_info_row_2_col[Y].info[6][row_id%ROW_LENGTH] = info0_6;
	d_info_row_2_col[Y].info[7][row_id%ROW_LENGTH] = info0_7;
	d_info_row_2_col[Y].info[8][row_id%ROW_LENGTH] = info0_8;
	d_info_row_2_col[Y].info[9][row_id%ROW_LENGTH] = info0_9;




	/**************************ROW_2********************************/
	row_id = row_offset - I * 4 * 768 + 2 * 256 + tid;
	info_sum0 = 0; symbol0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 3 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 5 + (97 + tid) & 255];
		info0_2 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 6 + (91 + tid) & 255];
		info0_3 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 0 + (17 + tid) & 255];
		info0_4 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 4 + (0 + tid) & 255];
		info0_5 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 6 + (164 + tid) & 255];
		info0_6 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 0 + (11 + tid) & 255];
		info0_7 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 6 + (125 + tid) & 255];
		info0_8 = buf_d_channel_info[Y][256 * 2 + (237 + tid) & 255];
		info0_9 = buf_d_channel_info[Y][256 * 6 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 256 + ((d_matrix_node_c[block_id].order[3].number + tid) & 255)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 256 + ((d_matrix_node_c[block_id].order[4].number + tid) & 255)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 256 + ((d_matrix_node_c[block_id].order[5].number + tid) & 255)) % COL_LENGTH;
		node_pos0_6 = (offset + d_matrix_node_c[block_id].order[6].col * 256 + ((d_matrix_node_c[block_id].order[6].number + tid) & 255)) % COL_LENGTH;
		node_pos0_7 = (offset + d_matrix_node_c[block_id].order[7].col * 256 + ((d_matrix_node_c[block_id].order[7].number + tid) & 255)) % COL_LENGTH;
		node_pos0_8 = (offset + d_matrix_node_c[block_id].order[8].col * 256 + ((d_matrix_node_c[block_id].order[8].number + tid) & 255)) % COL_LENGTH;
		node_pos0_9 = (offset + d_matrix_node_c[block_id].order[9].col * 256 + ((d_matrix_node_c[block_id].order[9].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;
		node_cpos6 = (d_matrix_node_c[block_id].order[6].col_2_row) % COL_LENGTH;
		node_cpos7 = (d_matrix_node_c[block_id].order[7].col_2_row) % COL_LENGTH;
		node_cpos8 = (d_matrix_node_c[block_id].order[8].col_2_row) % COL_LENGTH;
		node_cpos9 = (d_matrix_node_c[block_id].order[9].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];
		info0_6 = d_info_col_2_row[Y].info[node_cpos6][node_pos0_6];
		info0_7 = d_info_col_2_row[Y].info[node_cpos7][node_pos0_7];
		info0_8 = d_info_col_2_row[Y].info[node_cpos8][node_pos0_8];
		info0_9 = d_info_col_2_row[Y].info[node_cpos9][node_pos0_9];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info_symbol6 = ((unsigned int)(int)(info0_6 - 1)) >> 31;
	info0_6 = dLtanh(info0_6);
	info_sum_symbol0 ^= info_symbol6;
	info_sum0 += info0_6;

	info_symbol7 = ((unsigned int)(int)(info0_7 - 1)) >> 31;
	info0_7 = dLtanh(info0_7);
	info_sum_symbol0 ^= info_symbol7;
	info_sum0 += info0_7;

	info_symbol8 = ((unsigned int)(int)(info0_8 - 1)) >> 31;
	info0_8 = dLtanh(info0_8);
	info_sum_symbol0 ^= info_symbol8;
	info_sum0 += info0_8;

	info_symbol9 = ((unsigned int)(int)(info0_9 - 1)) >> 31;
	info0_9 = dLtanh(info0_9);
	info_sum_symbol0 ^= info_symbol9;
	info_sum0 += info0_9;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol5));

	info0_6 = info_sum0 - info0_6;
	info0_6 = dLtanh(info0_6);
	info0_6 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol6));

	info0_7 = info_sum0 - info0_7;
	info0_7 = dLtanh(info0_7);
	info0_7 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol7));

	info0_8 = info_sum0 - info0_8;
	info0_8 = dLtanh(info0_8);
	info0_8 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol8));

	info0_9 = info_sum0 - info0_9;
	info0_9 = dLtanh(info0_9);
	info0_9 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol9));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;
	d_info_row_2_col[Y].info[6][row_id%ROW_LENGTH] = info0_6;
	d_info_row_2_col[Y].info[7][row_id%ROW_LENGTH] = info0_7;
	d_info_row_2_col[Y].info[8][row_id%ROW_LENGTH] = info0_8;
	d_info_row_2_col[Y].info[9][row_id%ROW_LENGTH] = info0_9;


	__syncthreads();




	/***************************************************************/
	/***********************COL_UPDATE******************************/
	/***************************************************************/

	/**************************COL_0********************************/
	col_id = col_offset - I * 4 * 1792 + 0 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;//256;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][0 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(0 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
	}



	/**************************COL_1********************************/
	col_id = col_offset - I * 4 * 1792 + 1 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][1 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(1 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
	}



	/**************************COL_2********************************/
	col_id = col_offset - I * 4 * 1792 + 2 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][2 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(2 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
	}



	/**************************COL_3********************************/
	col_id = col_offset - I * 4 * 1792 + 3 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][3 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(3 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
	}



	/**************************COL_4********************************/
	col_id = col_offset - I * 4 * 1792 + 4 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][4 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(4 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
	}



	/**************************COL_5********************************/
	col_id = col_offset - I * 4 * 1792 + 5 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][5 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(5 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
	}



	/**************************COL_6********************************/
	col_id = col_offset - I * 4 * 1792 + 6 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][6 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(6 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
	}



	__syncthreads();


};

__global__ void dUpdate2(__int64 time_count, info_ch *d_channel_info, int *d_decoded_word, INFO_COL *d_info_col_2_row, INFO_ROW *d_info_row_2_col
	, buf_info_ch *buf_d_channel_info, BUF_INFO_COL *buf_d_info_col_2_row, int STREAM_COUNT
	)
{
	int tid = threadIdx.x;
	int I = blockIdx.x;
	int Y = blockIdx.y;

	int row_offset = time_count * 768;
	int col_offset = (time_count - 3) * 1792;
	int row_id = 0;
	int col_id = 0;
	register int node_pos0_0, node_pos0_1, node_pos0_2, node_pos0_3, node_pos0_4, node_pos0_5, node_pos0_6, node_pos0_7, node_pos0_8, node_pos0_9;
	register int node_cpos0, node_cpos1, node_cpos2, node_cpos3, node_cpos4, node_cpos5, node_cpos6, node_cpos7, node_cpos8, node_cpos9;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	register float info0_0, info0_1, info0_2, info0_3, info0_4, info0_5, info0_6, info0_7, info0_8, info0_9; //info = channel col to row info
	register int info_symbol0, info_symbol1, info_symbol2, info_symbol3, info_symbol4, info_symbol5, info_symbol6, info_symbol7, info_symbol8, info_symbol9; //positive = 0, negative = 1
	register float info_sum0;
	register int info_sum_symbol0;
	register int symbol0;
	register float info_channel0;
	int node_pos_2 = 0;


	/***************************************************************/
	/***********************ROW_UPDATE******************************/
	/***************************************************************/

	/**************************ROW_0********************************/
	row_id = row_offset - I * 4 * 768 + 0 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 2 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[Y][256 * 3 + (0 + tid) & 255];
		info0_2 = buf_d_channel_info[Y][256 * 6 + (160 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;



	/**************************ROW_1********************************/
	row_id = row_offset - I * 4 * 768 + 1 * 256 + tid;
	info_sum0 = 0; symbol0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 3 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[Y][256 * 3 + (0 + tid) & 255];
		info0_2 = buf_d_channel_info[Y][256 * 4 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;




	/**************************ROW_2********************************/
	row_id = row_offset - I * 4 * 768 + 2 * 256 + tid;
	info_sum0 = 0; symbol0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 4 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[Y][256 * 4 + (0 + tid) & 255];
		info0_2 = buf_d_channel_info[Y][256 * 5 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;


	__syncthreads();




	/***************************************************************/
	/***********************COL_UPDATE******************************/
	/***************************************************************/

	/**************************COL_0********************************/
	col_id = col_offset - I * 4 * 1792 + 0 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][0 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(0 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
	}



	/**************************COL_1********************************/
	col_id = col_offset - I * 4 * 1792 + 1 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][1 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];


	info_sum0 = info0_0 + info0_1;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(1 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;

	}



	/**************************COL_2********************************/
	col_id = col_offset - I * 4 * 1792 + 2 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][2 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];


	info_sum0 = info0_0 + info0_1;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(2 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;

	}



	/**************************COL_3********************************/
	col_id = col_offset - I * 4 * 1792 + 3 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][3 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];


	info_sum0 = info0_0 + info0_1;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(3 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;

	}



	/**************************COL_4********************************/
	col_id = col_offset - I * 4 * 1792 + 4 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][4 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];


	info_sum0 = info0_0 + info0_1;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(4 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;

	}



	/**************************COL_5********************************/
	col_id = col_offset - I * 4 * 1792 + 5 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][5 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];


	info_sum0 = info0_0 + info0_1 + info0_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(5 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

	}



	/**************************COL_6********************************/
	col_id = col_offset - I * 4 * 1792 + 6 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][6 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];


	info_sum0 = info0_0 + info0_1 + info0_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(6 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

	}



	__syncthreads();


};


__global__ void dUpdate3(__int64 time_count, info_ch *d_channel_info, int *d_decoded_word, INFO_COL *d_info_col_2_row, INFO_ROW *d_info_row_2_col
	, buf_info_ch *buf_d_channel_info, BUF_INFO_COL *buf_d_info_col_2_row, int STREAM_COUNT
	)
{
	int tid = threadIdx.x;
	int I = blockIdx.x;
	int Y = blockIdx.y;

	int row_offset = time_count * 768;
	int col_offset = (time_count - 3) * 1792;
	int row_id = 0;
	int col_id = 0;
	register int node_pos0_0, node_pos0_1, node_pos0_2, node_pos0_3, node_pos0_4, node_pos0_5, node_pos0_6, node_pos0_7, node_pos0_8, node_pos0_9;
	register int node_cpos0, node_cpos1, node_cpos2, node_cpos3, node_cpos4, node_cpos5, node_cpos6, node_cpos7, node_cpos8, node_cpos9;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	register float info0_0, info0_1, info0_2, info0_3, info0_4, info0_5, info0_6, info0_7, info0_8, info0_9;
	register int info_symbol0, info_symbol1, info_symbol2, info_symbol3, info_symbol4, info_symbol5, info_symbol6, info_symbol7, info_symbol8, info_symbol9; //positive = 0, negative = 1
	register float info_sum0;
	register int info_sum_symbol0;
	register int symbol0;
	register float info_channel0;
	int node_pos_2 = 0;


	/***************************************************************/
	/***********************ROW_UPDATE******************************/
	/***************************************************************/

	/**************************ROW_0********************************/
	row_id = row_offset - I * 4 * 768 + 0 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 5 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[Y][256 * 5 + (0 + tid) & 255];
		info0_2 = buf_d_channel_info[Y][256 * 6 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;



	/**************************ROW_1********************************/
	row_id = row_offset - I * 4 * 768 + 1 * 256 + tid;
	info_sum0 = 0; symbol0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 6 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 3 + (241 + tid) & 255];
		info0_2 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 4 + (185 + tid) & 255];
		info0_3 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 5 + (251 + tid) & 255];
		info0_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 0 + (248 + tid) & 255];
		info0_5 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 1 + (12 + tid) & 255];
		info0_6 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 2 + (11 + tid) & 255];
		info0_7 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 4 + (0 + tid) & 255];
		info0_8 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 1 + (0 + tid) & 255];
		info0_9 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 5 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 256 + ((d_matrix_node_c[block_id].order[3].number + tid) & 255)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 256 + ((d_matrix_node_c[block_id].order[4].number + tid) & 255)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 256 + ((d_matrix_node_c[block_id].order[5].number + tid) & 255)) % COL_LENGTH;
		node_pos0_6 = (offset + d_matrix_node_c[block_id].order[6].col * 256 + ((d_matrix_node_c[block_id].order[6].number + tid) & 255)) % COL_LENGTH;
		node_pos0_7 = (offset + d_matrix_node_c[block_id].order[7].col * 256 + ((d_matrix_node_c[block_id].order[7].number + tid) & 255)) % COL_LENGTH;
		node_pos0_8 = (offset + d_matrix_node_c[block_id].order[8].col * 256 + ((d_matrix_node_c[block_id].order[8].number + tid) & 255)) % COL_LENGTH;
		node_pos0_9 = (offset + d_matrix_node_c[block_id].order[9].col * 256 + ((d_matrix_node_c[block_id].order[9].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;
		node_cpos6 = (d_matrix_node_c[block_id].order[6].col_2_row) % COL_LENGTH;
		node_cpos7 = (d_matrix_node_c[block_id].order[7].col_2_row) % COL_LENGTH;
		node_cpos8 = (d_matrix_node_c[block_id].order[8].col_2_row) % COL_LENGTH;
		node_cpos9 = (d_matrix_node_c[block_id].order[9].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];
		info0_6 = d_info_col_2_row[Y].info[node_cpos6][node_pos0_6];
		info0_7 = d_info_col_2_row[Y].info[node_cpos7][node_pos0_7];
		info0_8 = d_info_col_2_row[Y].info[node_cpos8][node_pos0_8];
		info0_9 = d_info_col_2_row[Y].info[node_cpos9][node_pos0_9];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info_symbol6 = ((unsigned int)(int)(info0_6 - 1)) >> 31;
	info0_6 = dLtanh(info0_6);
	info_sum_symbol0 ^= info_symbol6;
	info_sum0 += info0_6;

	info_symbol7 = ((unsigned int)(int)(info0_7 - 1)) >> 31;
	info0_7 = dLtanh(info0_7);
	info_sum_symbol0 ^= info_symbol7;
	info_sum0 += info0_7;

	info_symbol8 = ((unsigned int)(int)(info0_8 - 1)) >> 31;
	info0_8 = dLtanh(info0_8);
	info_sum_symbol0 ^= info_symbol8;
	info_sum0 += info0_8;

	info_symbol9 = ((unsigned int)(int)(info0_9 - 1)) >> 31;
	info0_9 = dLtanh(info0_9);
	info_sum_symbol0 ^= info_symbol9;
	info_sum0 += info0_9;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol5));

	info0_6 = info_sum0 - info0_6;
	info0_6 = dLtanh(info0_6);
	info0_6 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol6));

	info0_7 = info_sum0 - info0_7;
	info0_7 = dLtanh(info0_7);
	info0_7 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol7));

	info0_8 = info_sum0 - info0_8;
	info0_8 = dLtanh(info0_8);
	info0_8 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol8));

	info0_9 = info_sum0 - info0_9;
	info0_9 = dLtanh(info0_9);
	info0_9 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol9));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;
	d_info_row_2_col[Y].info[6][row_id%ROW_LENGTH] = info0_6;
	d_info_row_2_col[Y].info[7][row_id%ROW_LENGTH] = info0_7;
	d_info_row_2_col[Y].info[8][row_id%ROW_LENGTH] = info0_8;
	d_info_row_2_col[Y].info[9][row_id%ROW_LENGTH] = info0_9;




	/**************************ROW_2********************************/
	row_id = row_offset - I * 4 * 768 + 2 * 256 + tid;
	info_sum0 = 0; symbol0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 0 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 4 + (182 + tid) & 255];
		info0_2 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 5 + (249 + tid) & 255];
		info0_3 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 6 + (65 + tid) & 255];
		info0_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 1 + (55 + tid) & 255];
		info0_5 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 2 + (12 + tid) & 255];
		info0_6 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 3 + (227 + tid) & 255];
		info0_7 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 5 + (0 + tid) & 255];
		info0_8 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 2 + (0 + tid) & 255];
		info0_9 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 6 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 256 + ((d_matrix_node_c[block_id].order[3].number + tid) & 255)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 256 + ((d_matrix_node_c[block_id].order[4].number + tid) & 255)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 256 + ((d_matrix_node_c[block_id].order[5].number + tid) & 255)) % COL_LENGTH;
		node_pos0_6 = (offset + d_matrix_node_c[block_id].order[6].col * 256 + ((d_matrix_node_c[block_id].order[6].number + tid) & 255)) % COL_LENGTH;
		node_pos0_7 = (offset + d_matrix_node_c[block_id].order[7].col * 256 + ((d_matrix_node_c[block_id].order[7].number + tid) & 255)) % COL_LENGTH;
		node_pos0_8 = (offset + d_matrix_node_c[block_id].order[8].col * 256 + ((d_matrix_node_c[block_id].order[8].number + tid) & 255)) % COL_LENGTH;
		node_pos0_9 = (offset + d_matrix_node_c[block_id].order[9].col * 256 + ((d_matrix_node_c[block_id].order[9].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;
		node_cpos6 = (d_matrix_node_c[block_id].order[6].col_2_row) % COL_LENGTH;
		node_cpos7 = (d_matrix_node_c[block_id].order[7].col_2_row) % COL_LENGTH;
		node_cpos8 = (d_matrix_node_c[block_id].order[8].col_2_row) % COL_LENGTH;
		node_cpos9 = (d_matrix_node_c[block_id].order[9].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];
		info0_6 = d_info_col_2_row[Y].info[node_cpos6][node_pos0_6];
		info0_7 = d_info_col_2_row[Y].info[node_cpos7][node_pos0_7];
		info0_8 = d_info_col_2_row[Y].info[node_cpos8][node_pos0_8];
		info0_9 = d_info_col_2_row[Y].info[node_cpos9][node_pos0_9];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info_symbol6 = ((unsigned int)(int)(info0_6 - 1)) >> 31;
	info0_6 = dLtanh(info0_6);
	info_sum_symbol0 ^= info_symbol6;
	info_sum0 += info0_6;

	info_symbol7 = ((unsigned int)(int)(info0_7 - 1)) >> 31;
	info0_7 = dLtanh(info0_7);
	info_sum_symbol0 ^= info_symbol7;
	info_sum0 += info0_7;

	info_symbol8 = ((unsigned int)(int)(info0_8 - 1)) >> 31;
	info0_8 = dLtanh(info0_8);
	info_sum_symbol0 ^= info_symbol8;
	info_sum0 += info0_8;

	info_symbol9 = ((unsigned int)(int)(info0_9 - 1)) >> 31;
	info0_9 = dLtanh(info0_9);
	info_sum_symbol0 ^= info_symbol9;
	info_sum0 += info0_9;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol5));

	info0_6 = info_sum0 - info0_6;
	info0_6 = dLtanh(info0_6);
	info0_6 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol6));

	info0_7 = info_sum0 - info0_7;
	info0_7 = dLtanh(info0_7);
	info0_7 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol7));

	info0_8 = info_sum0 - info0_8;
	info0_8 = dLtanh(info0_8);
	info0_8 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol8));

	info0_9 = info_sum0 - info0_9;
	info0_9 = dLtanh(info0_9);
	info0_9 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol9));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;
	d_info_row_2_col[Y].info[6][row_id%ROW_LENGTH] = info0_6;
	d_info_row_2_col[Y].info[7][row_id%ROW_LENGTH] = info0_7;
	d_info_row_2_col[Y].info[8][row_id%ROW_LENGTH] = info0_8;
	d_info_row_2_col[Y].info[9][row_id%ROW_LENGTH] = info0_9;


	__syncthreads();




	/***************************************************************/
	/***********************COL_UPDATE******************************/
	/***************************************************************/

	/**************************COL_0********************************/
	col_id = col_offset - I * 4 * 1792 + 0 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][0 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];


	info_sum0 = info0_0 + info0_1 + info0_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(0 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

	}



	/**************************COL_1********************************/
	col_id = col_offset - I * 4 * 1792 + 1 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][1 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];


	info_sum0 = info0_0 + info0_1 + info0_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(1 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

	}



	/**************************COL_2********************************/
	col_id = col_offset - I * 4 * 1792 + 2 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][2 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];


	info_sum0 = info0_0;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(2 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;

	}



	/**************************COL_3********************************/
	col_id = col_offset - I * 4 * 1792 + 3 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][3 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];


	info_sum0 = info0_0;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(3 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;

	}



	/**************************COL_4********************************/
	col_id = col_offset - I * 4 * 1792 + 4 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][4 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];


	info_sum0 = info0_0;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(4 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;

	}



	/**************************COL_5********************************/
	col_id = col_offset - I * 4 * 1792 + 5 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][5 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];


	info_sum0 = info0_0;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(5 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;

	}



	/**************************COL_6********************************/
	col_id = col_offset - I * 4 * 1792 + 6 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][6 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];


	info_sum0 = info0_0 + info0_1 + info0_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(6 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

	}



	__syncthreads();


};


__global__ void dUpdate4(__int64 time_count, info_ch *d_channel_info, int *d_decoded_word, INFO_COL *d_info_col_2_row, INFO_ROW *d_info_row_2_col
	, buf_info_ch *buf_d_channel_info, BUF_INFO_COL *buf_d_info_col_2_row, int STREAM_COUNT
	)
{
	int tid = threadIdx.x;
	int I = blockIdx.x;
	int Y = blockIdx.y;

	int row_offset = time_count * 768;
	int col_offset = (time_count - 3) * 1792;
	int row_id = 0;
	int col_id = 0;
	register int node_pos0_0, node_pos0_1, node_pos0_2, node_pos0_3, node_pos0_4, node_pos0_5, node_pos0_6, node_pos0_7, node_pos0_8, node_pos0_9;
	register int node_cpos0, node_cpos1, node_cpos2, node_cpos3, node_cpos4, node_cpos5, node_cpos6, node_cpos7, node_cpos8, node_cpos9;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	register float info0_0, info0_1, info0_2, info0_3, info0_4, info0_5, info0_6, info0_7, info0_8, info0_9; //info = channel col to row info
	register int info_symbol0, info_symbol1, info_symbol2, info_symbol3, info_symbol4, info_symbol5, info_symbol6, info_symbol7, info_symbol8, info_symbol9; //positive = 0, negative = 1
	register float info_sum0;
	register int info_sum_symbol0;
	register int symbol0;
	register float info_channel0;
	int node_pos_2 = 0;


	/***************************************************************/
	/***********************ROW_UPDATE******************************/
	/***************************************************************/

	/**************************ROW_0********************************/
	row_id = row_offset - I * 4 * 768 + 0 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 1 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 3 + (214 + tid) & 255];
		info0_2 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 5 + (35 + tid) & 255];
		info0_3 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 6 + (167 + tid) & 255];
		info0_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 0 + (23 + tid) & 255];
		info0_5 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 2 + (147 + tid) & 255];
		info0_6 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 3 + (54 + tid) & 255];
		info0_7 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 6 + (0 + tid) & 255];
		info0_8 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 3 + (0 + tid) & 255];
		info0_9 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 0 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 256 + ((d_matrix_node_c[block_id].order[3].number + tid) & 255)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 256 + ((d_matrix_node_c[block_id].order[4].number + tid) & 255)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 256 + ((d_matrix_node_c[block_id].order[5].number + tid) & 255)) % COL_LENGTH;
		node_pos0_6 = (offset + d_matrix_node_c[block_id].order[6].col * 256 + ((d_matrix_node_c[block_id].order[6].number + tid) & 255)) % COL_LENGTH;
		node_pos0_7 = (offset + d_matrix_node_c[block_id].order[7].col * 256 + ((d_matrix_node_c[block_id].order[7].number + tid) & 255)) % COL_LENGTH;
		node_pos0_8 = (offset + d_matrix_node_c[block_id].order[8].col * 256 + ((d_matrix_node_c[block_id].order[8].number + tid) & 255)) % COL_LENGTH;
		node_pos0_9 = (offset + d_matrix_node_c[block_id].order[9].col * 256 + ((d_matrix_node_c[block_id].order[9].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;
		node_cpos6 = (d_matrix_node_c[block_id].order[6].col_2_row) % COL_LENGTH;
		node_cpos7 = (d_matrix_node_c[block_id].order[7].col_2_row) % COL_LENGTH;
		node_cpos8 = (d_matrix_node_c[block_id].order[8].col_2_row) % COL_LENGTH;
		node_cpos9 = (d_matrix_node_c[block_id].order[9].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];
		info0_6 = d_info_col_2_row[Y].info[node_cpos6][node_pos0_6];
		info0_7 = d_info_col_2_row[Y].info[node_cpos7][node_pos0_7];
		info0_8 = d_info_col_2_row[Y].info[node_cpos8][node_pos0_8];
		info0_9 = d_info_col_2_row[Y].info[node_cpos9][node_pos0_9];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info_symbol6 = ((unsigned int)(int)(info0_6 - 1)) >> 31;
	info0_6 = dLtanh(info0_6);
	info_sum_symbol0 ^= info_symbol6;
	info_sum0 += info0_6;

	info_symbol7 = ((unsigned int)(int)(info0_7 - 1)) >> 31;
	info0_7 = dLtanh(info0_7);
	info_sum_symbol0 ^= info_symbol7;
	info_sum0 += info0_7;

	info_symbol8 = ((unsigned int)(int)(info0_8 - 1)) >> 31;
	info0_8 = dLtanh(info0_8);
	info_sum_symbol0 ^= info_symbol8;
	info_sum0 += info0_8;

	info_symbol9 = ((unsigned int)(int)(info0_9 - 1)) >> 31;
	info0_9 = dLtanh(info0_9);
	info_sum_symbol0 ^= info_symbol9;
	info_sum0 += info0_9;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol5));

	info0_6 = info_sum0 - info0_6;
	info0_6 = dLtanh(info0_6);
	info0_6 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol6));

	info0_7 = info_sum0 - info0_7;
	info0_7 = dLtanh(info0_7);
	info0_7 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol7));

	info0_8 = info_sum0 - info0_8;
	info0_8 = dLtanh(info0_8);
	info0_8 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol8));

	info0_9 = info_sum0 - info0_9;
	info0_9 = dLtanh(info0_9);
	info0_9 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol9));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;
	d_info_row_2_col[Y].info[6][row_id%ROW_LENGTH] = info0_6;
	d_info_row_2_col[Y].info[7][row_id%ROW_LENGTH] = info0_7;
	d_info_row_2_col[Y].info[8][row_id%ROW_LENGTH] = info0_8;
	d_info_row_2_col[Y].info[9][row_id%ROW_LENGTH] = info0_9;



	/**************************ROW_1********************************/
	row_id = row_offset - I * 4 * 768 + 1 * 256 + tid;
	info_sum0 = 0; symbol0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 2 + (0 + tid) & 255];
		info0_1 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 3 + (7 + tid) & 255];
		info0_2 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 4 + (31 + tid) & 255];
		info0_3 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 6 + (162 + tid) & 255];
		info0_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 0 + (99 + tid) & 255];
		info0_5 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 1 + (105 + tid) & 255];
		info0_6 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 3 + (133 + tid) & 255];
		info0_7 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 0 + (0 + tid) & 255];
		info0_8 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 4 + (0 + tid) & 255];
		info0_9 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 1 + (0 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 256 + ((d_matrix_node_c[block_id].order[3].number + tid) & 255)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 256 + ((d_matrix_node_c[block_id].order[4].number + tid) & 255)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 256 + ((d_matrix_node_c[block_id].order[5].number + tid) & 255)) % COL_LENGTH;
		node_pos0_6 = (offset + d_matrix_node_c[block_id].order[6].col * 256 + ((d_matrix_node_c[block_id].order[6].number + tid) & 255)) % COL_LENGTH;
		node_pos0_7 = (offset + d_matrix_node_c[block_id].order[7].col * 256 + ((d_matrix_node_c[block_id].order[7].number + tid) & 255)) % COL_LENGTH;
		node_pos0_8 = (offset + d_matrix_node_c[block_id].order[8].col * 256 + ((d_matrix_node_c[block_id].order[8].number + tid) & 255)) % COL_LENGTH;
		node_pos0_9 = (offset + d_matrix_node_c[block_id].order[9].col * 256 + ((d_matrix_node_c[block_id].order[9].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;
		node_cpos6 = (d_matrix_node_c[block_id].order[6].col_2_row) % COL_LENGTH;
		node_cpos7 = (d_matrix_node_c[block_id].order[7].col_2_row) % COL_LENGTH;
		node_cpos8 = (d_matrix_node_c[block_id].order[8].col_2_row) % COL_LENGTH;
		node_cpos9 = (d_matrix_node_c[block_id].order[9].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];
		info0_6 = d_info_col_2_row[Y].info[node_cpos6][node_pos0_6];
		info0_7 = d_info_col_2_row[Y].info[node_cpos7][node_pos0_7];
		info0_8 = d_info_col_2_row[Y].info[node_cpos8][node_pos0_8];
		info0_9 = d_info_col_2_row[Y].info[node_cpos9][node_pos0_9];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info_symbol6 = ((unsigned int)(int)(info0_6 - 1)) >> 31;
	info0_6 = dLtanh(info0_6);
	info_sum_symbol0 ^= info_symbol6;
	info_sum0 += info0_6;

	info_symbol7 = ((unsigned int)(int)(info0_7 - 1)) >> 31;
	info0_7 = dLtanh(info0_7);
	info_sum_symbol0 ^= info_symbol7;
	info_sum0 += info0_7;

	info_symbol8 = ((unsigned int)(int)(info0_8 - 1)) >> 31;
	info0_8 = dLtanh(info0_8);
	info_sum_symbol0 ^= info_symbol8;
	info_sum0 += info0_8;

	info_symbol9 = ((unsigned int)(int)(info0_9 - 1)) >> 31;
	info0_9 = dLtanh(info0_9);
	info_sum_symbol0 ^= info_symbol9;
	info_sum0 += info0_9;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol5));

	info0_6 = info_sum0 - info0_6;
	info0_6 = dLtanh(info0_6);
	info0_6 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol6));

	info0_7 = info_sum0 - info0_7;
	info0_7 = dLtanh(info0_7);
	info0_7 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol7));

	info0_8 = info_sum0 - info0_8;
	info0_8 = dLtanh(info0_8);
	info0_8 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol8));

	info0_9 = info_sum0 - info0_9;
	info0_9 = dLtanh(info0_9);
	info0_9 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol9));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;
	d_info_row_2_col[Y].info[6][row_id%ROW_LENGTH] = info0_6;
	d_info_row_2_col[Y].info[7][row_id%ROW_LENGTH] = info0_7;
	d_info_row_2_col[Y].info[8][row_id%ROW_LENGTH] = info0_8;
	d_info_row_2_col[Y].info[9][row_id%ROW_LENGTH] = info0_9;




	/**************************ROW_2********************************/
	row_id = row_offset - I * 4 * 768 + 2 * 256 + tid;
	info_sum0 = 0; symbol0 = 0;
	info_sum_symbol0 = 0;

	block_id = (row_id % 3072) >> 8;
	offset = ((row_id / 768) - ROW_OFFSET) * 1792;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 0 + (184 + tid) & 255];
		info0_1 = buf_d_channel_info[0 * STREAM_COUNT + Y][256 * 3 + (0 + tid) & 255];
		info0_2 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 0 + (0 + tid) & 255];
		info0_3 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 4 + (66 + tid) & 255];
		info0_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][256 * 6 + (173 + tid) & 255];
		info0_5 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 0 + (42 + tid) & 255];
		info0_6 = buf_d_channel_info[2 * STREAM_COUNT + Y][256 * 1 + (0 + tid) & 255];
		info0_7 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 0 + (209 + tid) & 255];
		info0_8 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 1 + (103 + tid) & 255];
		info0_9 = buf_d_channel_info[3 * STREAM_COUNT + Y][256 * 6 + (90 + tid) & 255];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 256 + ((d_matrix_node_c[block_id].order[0].number + tid) & 255)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 256 + ((d_matrix_node_c[block_id].order[1].number + tid) & 255)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 256 + ((d_matrix_node_c[block_id].order[2].number + tid) & 255)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 256 + ((d_matrix_node_c[block_id].order[3].number + tid) & 255)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 256 + ((d_matrix_node_c[block_id].order[4].number + tid) & 255)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 256 + ((d_matrix_node_c[block_id].order[5].number + tid) & 255)) % COL_LENGTH;
		node_pos0_6 = (offset + d_matrix_node_c[block_id].order[6].col * 256 + ((d_matrix_node_c[block_id].order[6].number + tid) & 255)) % COL_LENGTH;
		node_pos0_7 = (offset + d_matrix_node_c[block_id].order[7].col * 256 + ((d_matrix_node_c[block_id].order[7].number + tid) & 255)) % COL_LENGTH;
		node_pos0_8 = (offset + d_matrix_node_c[block_id].order[8].col * 256 + ((d_matrix_node_c[block_id].order[8].number + tid) & 255)) % COL_LENGTH;
		node_pos0_9 = (offset + d_matrix_node_c[block_id].order[9].col * 256 + ((d_matrix_node_c[block_id].order[9].number + tid) & 255)) % COL_LENGTH;


		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;
		node_cpos6 = (d_matrix_node_c[block_id].order[6].col_2_row) % COL_LENGTH;
		node_cpos7 = (d_matrix_node_c[block_id].order[7].col_2_row) % COL_LENGTH;
		node_cpos8 = (d_matrix_node_c[block_id].order[8].col_2_row) % COL_LENGTH;
		node_cpos9 = (d_matrix_node_c[block_id].order[9].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];
		info0_6 = d_info_col_2_row[Y].info[node_cpos6][node_pos0_6];
		info0_7 = d_info_col_2_row[Y].info[node_cpos7][node_pos0_7];
		info0_8 = d_info_col_2_row[Y].info[node_cpos8][node_pos0_8];
		info0_9 = d_info_col_2_row[Y].info[node_cpos9][node_pos0_9];
	}

	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info_symbol6 = ((unsigned int)(int)(info0_6 - 1)) >> 31;
	info0_6 = dLtanh(info0_6);
	info_sum_symbol0 ^= info_symbol6;
	info_sum0 += info0_6;

	info_symbol7 = ((unsigned int)(int)(info0_7 - 1)) >> 31;
	info0_7 = dLtanh(info0_7);
	info_sum_symbol0 ^= info_symbol7;
	info_sum0 += info0_7;

	info_symbol8 = ((unsigned int)(int)(info0_8 - 1)) >> 31;
	info0_8 = dLtanh(info0_8);
	info_sum_symbol0 ^= info_symbol8;
	info_sum0 += info0_8;

	info_symbol9 = ((unsigned int)(int)(info0_9 - 1)) >> 31;
	info0_9 = dLtanh(info0_9);
	info_sum_symbol0 ^= info_symbol9;
	info_sum0 += info0_9;


	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol5));

	info0_6 = info_sum0 - info0_6;
	info0_6 = dLtanh(info0_6);
	info0_6 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol6));

	info0_7 = info_sum0 - info0_7;
	info0_7 = dLtanh(info0_7);
	info0_7 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol7));

	info0_8 = info_sum0 - info0_8;
	info0_8 = dLtanh(info0_8);
	info0_8 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol8));

	info0_9 = info_sum0 - info0_9;
	info0_9 = dLtanh(info0_9);
	info0_9 *= (1 - 2 * (info_sum_symbol0 ^ info_symbol9));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;
	d_info_row_2_col[Y].info[6][row_id%ROW_LENGTH] = info0_6;
	d_info_row_2_col[Y].info[7][row_id%ROW_LENGTH] = info0_7;
	d_info_row_2_col[Y].info[8][row_id%ROW_LENGTH] = info0_8;
	d_info_row_2_col[Y].info[9][row_id%ROW_LENGTH] = info0_9;


	__syncthreads();




	/***************************************************************/
	/***********************COL_UPDATE******************************/
	/***************************************************************/

	/**************************COL_0********************************/
	col_id = col_offset - I * 4 * 1792 + 0 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[0 * STREAM_COUNT + Y][0 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];


	info_sum0 = info0_0 + info0_1 + info0_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(0 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

	}



	/**************************COL_1********************************/
	col_id = col_offset - I * 4 * 1792 + 1 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[0 * STREAM_COUNT + Y][1 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];


	info_sum0 = info0_0 + info0_1 + info0_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(1 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

	}



	/**************************COL_2********************************/
	col_id = col_offset - I * 4 * 1792 + 2 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[0 * STREAM_COUNT + Y][2 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];


	info_sum0 = info0_0 + info0_1 + info0_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(2 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

	}



	/**************************COL_3********************************/
	col_id = col_offset - I * 4 * 1792 + 3 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[0 * STREAM_COUNT + Y][3 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_4 = (offset + d_matrix_node_v[block_id].order[4].row * 256 + ((tid - d_matrix_node_c[block_id].order[4].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_5 = (offset + d_matrix_node_v[block_id].order[5].row * 256 + ((tid - d_matrix_node_c[block_id].order[5].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;
	node_cpos4 = (d_matrix_node_v[block_id].order[4].row_2_col) % ROW_LENGTH;
	node_cpos5 = (d_matrix_node_v[block_id].order[5].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];
	info0_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos0_4];
	info0_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos0_5];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3 + info0_4 + info0_5;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(3 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;
		info0_4 = info_sum0 - info0_4 + info_channel0;
		info0_5 = info_sum0 - info0_5 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
		d_info_col_2_row[Y].info[4][col_id%COL_LENGTH] = info0_4;
		d_info_col_2_row[Y].info[5][col_id%COL_LENGTH] = info0_5;

	}



	/**************************COL_4********************************/
	col_id = col_offset - I * 4 * 1792 + 4 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[0 * STREAM_COUNT + Y][4 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_4 = (offset + d_matrix_node_v[block_id].order[4].row * 256 + ((tid - d_matrix_node_c[block_id].order[4].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_5 = (offset + d_matrix_node_v[block_id].order[5].row * 256 + ((tid - d_matrix_node_c[block_id].order[5].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;
	node_cpos4 = (d_matrix_node_v[block_id].order[4].row_2_col) % ROW_LENGTH;
	node_cpos5 = (d_matrix_node_v[block_id].order[5].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];
	info0_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos0_4];
	info0_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos0_5];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3 + info0_4 + info0_5;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(4 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;
		info0_4 = info_sum0 - info0_4 + info_channel0;
		info0_5 = info_sum0 - info0_5 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
		d_info_col_2_row[Y].info[4][col_id%COL_LENGTH] = info0_4;
		d_info_col_2_row[Y].info[5][col_id%COL_LENGTH] = info0_5;

	}



	/**************************COL_5********************************/
	col_id = col_offset - I * 4 * 1792 + 5 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[0 * STREAM_COUNT + Y][5 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_4 = (offset + d_matrix_node_v[block_id].order[4].row * 256 + ((tid - d_matrix_node_c[block_id].order[4].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_5 = (offset + d_matrix_node_v[block_id].order[5].row * 256 + ((tid - d_matrix_node_c[block_id].order[5].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;
	node_cpos4 = (d_matrix_node_v[block_id].order[4].row_2_col) % ROW_LENGTH;
	node_cpos5 = (d_matrix_node_v[block_id].order[5].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];
	info0_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos0_4];
	info0_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos0_5];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3 + info0_4 + info0_5;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(5 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;
		info0_4 = info_sum0 - info0_4 + info_channel0;
		info0_5 = info_sum0 - info0_5 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
		d_info_col_2_row[Y].info[4][col_id%COL_LENGTH] = info0_4;
		d_info_col_2_row[Y].info[5][col_id%COL_LENGTH] = info0_5;

	}



	/**************************COL_6********************************/
	col_id = col_offset - I * 4 * 1792 + 6 * 256 + tid;
	info_sum0 = 0;
	info_sum_symbol0 = 0;
	info_channel0 = 0;

	block_id = (col_id % 7168) >> 8;
	offset = (col_id / 1792) * 768;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[0 * STREAM_COUNT + Y][6 * 256 + tid];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 256 + ((tid - d_matrix_node_c[block_id].order[0].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 256 + ((tid - d_matrix_node_c[block_id].order[1].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 256 + ((tid - d_matrix_node_c[block_id].order[2].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 256 + ((tid - d_matrix_node_c[block_id].order[3].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_4 = (offset + d_matrix_node_v[block_id].order[4].row * 256 + ((tid - d_matrix_node_c[block_id].order[4].number + 256) & 255)) % ROW_LENGTH;
	node_pos0_5 = (offset + d_matrix_node_v[block_id].order[5].row * 256 + ((tid - d_matrix_node_c[block_id].order[5].number + 256) & 255)) % ROW_LENGTH;


	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;
	node_cpos4 = (d_matrix_node_v[block_id].order[4].row_2_col) % ROW_LENGTH;
	node_cpos5 = (d_matrix_node_v[block_id].order[5].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];
	info0_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos0_4];
	info0_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos0_5];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3 + info0_4 + info0_5;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(6 * 256 + tid) % 1792 + Y * 1792] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;
		info0_4 = info_sum0 - info0_4 + info_channel0;
		info0_5 = info_sum0 - info0_5 + info_channel0;


		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
		d_info_col_2_row[Y].info[4][col_id%COL_LENGTH] = info0_4;
		d_info_col_2_row[Y].info[5][col_id%COL_LENGTH] = info0_5;

	}



	__syncthreads();


};


#endif
