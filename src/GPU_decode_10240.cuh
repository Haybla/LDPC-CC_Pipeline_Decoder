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

#ifdef CODE1

__noinline__ __device__ float dLtanh(float info)
{
	float y, x_abs, x_e;

	x_abs = (float)fabs(info);

	x_e = (float)__expf(-x_abs);
	y = (float)__logf((1 + x_e) / (1 - x_e));

	return (y);
};


__global__ void dUpdate1(int time_count, info_ch *d_channel_info, int *d_decoded_word, INFO_COL *d_info_col_2_row, INFO_ROW *d_info_row_2_col
	, buf_info_ch *buf_d_channel_info, BUF_INFO_COL *buf_d_info_col_2_row, int STREAM_COUNT
	)
{
	int tid = threadIdx.x;
	int I = blockIdx.x;
	int Y = blockIdx.y;

	int row_offset = time_count * 1536; 
	int col_offset = (time_count - 3) * 2560; 
	int row_id = 0;
	int col_id = 0;
	register int node_pos0_0, node_pos0_1, node_pos0_2, node_pos0_3, node_pos0_4, node_pos0_5;
	register int node_pos1_0, node_pos1_1, node_pos1_2, node_pos1_3, node_pos1_4, node_pos1_5;
	register int node_cpos0, node_cpos1, node_cpos2, node_cpos3, node_cpos4, node_cpos5;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	register float info0_0, info0_1, info0_2, info0_3, info0_4, info0_5; 
	register float info1_0, info1_1, info1_2, info1_3, info1_4, info1_5;
	register int info_symbol0, info_symbol1, info_symbol2, info_symbol3, info_symbol4, info_symbol5; 
	register float info_sum0, info_sum1;
	register int info_sum_symbol0, info_sum_symbol1;
	register int symbol0, symbol1;
	register float info_channel0, info_channel1;
	int node_pos_2 = 0;


	/***************************************************************/
	/***********************ROW_UPDATE******************************/
	/***************************************************************/

	/**************************ROW_0********************************/
	row_id = row_offset - I * 4 * 1536 + 0 * 512 + tid; 
	info_sum0 = 0; info_sum1 = 0; symbol0 = 0; symbol1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560; 


	if (I == 0)
	{
		info0_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 1 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 4 + (83 + tid) & 511];
		info0_2 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 2 + (260 + tid) & 511];
		info0_3 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 3 + (318 + tid) & 511];
		info0_4 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 4 + (382 + tid) & 511];
		info0_5 = buf_d_channel_info[Y][512 * 2 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 1 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 4 + (83 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 2 + (260 + tid + 256) & 511];
		info1_3 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 3 + (318 + tid + 256) & 511];
		info1_4 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 4 + (382 + tid + 256) & 511];
		info1_5 = buf_d_channel_info[Y][512 * 2 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid) & 511)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid) & 511)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
		info1_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos1_3];
		info1_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos1_4];
		info1_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos1_5];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	symbol0 |= info_symbol0 << 0;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	symbol0 |= info_symbol1 << 1;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	symbol0 |= info_symbol2 << 2;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	symbol0 |= info_symbol3 << 3;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	symbol0 |= info_symbol4 << 4;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	symbol0 |= info_symbol5 << 5;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 0 & 0x01)));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 1 & 0x01)));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 2 & 0x01)));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 3 & 0x01)));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 4 & 0x01)));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 5 & 0x01)));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	symbol1 |= info_symbol0 << 0;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	symbol1 |= info_symbol1 << 1;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	symbol1 |= info_symbol2 << 2;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;

	info_symbol3 = ((unsigned int)(int)(info1_3 - 1)) >> 31;
	symbol1 |= info_symbol3 << 3;
	info1_3 = dLtanh(info1_3);
	info_sum_symbol1 ^= info_symbol3;
	info_sum1 += info1_3;

	info_symbol4 = ((unsigned int)(int)(info1_4 - 1)) >> 31;
	symbol1 |= info_symbol4 << 4;
	info1_4 = dLtanh(info1_4);
	info_sum_symbol1 ^= info_symbol4;
	info_sum1 += info1_4;

	info_symbol5 = ((unsigned int)(int)(info1_5 - 1)) >> 31;
	symbol1 |= info_symbol5 << 5;
	info1_5 = dLtanh(info1_5);
	info_sum_symbol1 ^= info_symbol5;
	info_sum1 += info1_5;

	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 0 & 0x01)));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 1 & 0x01)));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 2 & 0x01)));

	info1_3 = info_sum1 - info1_3;
	info1_3 = dLtanh(info1_3);
	info1_3 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 3 & 0x01)));

	info1_4 = info_sum1 - info1_4;
	info1_4 = dLtanh(info1_4);
	info1_4 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 4 & 0x01)));

	info1_5 = info_sum1 - info1_5;
	info1_5 = dLtanh(info1_5);
	info1_5 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 5 & 0x01)));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;
	d_info_row_2_col[Y].info[3][(row_id + 256) % ROW_LENGTH] = info1_3;
	d_info_row_2_col[Y].info[4][(row_id + 256) % ROW_LENGTH] = info1_4;
	d_info_row_2_col[Y].info[5][(row_id + 256) % ROW_LENGTH] = info1_5;



	/**************************ROW_1********************************/
	row_id = row_offset - I * 4 * 1536 + 1 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0; symbol0 = 0; symbol1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 2 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 4 + (415 + tid) & 511];
		info0_2 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 0 + (403 + tid) & 511];
		info0_3 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 4 + (184 + tid) & 511];
		info0_4 = buf_d_channel_info[Y][512 * 0 + (279 + tid) & 511];
		info0_5 = buf_d_channel_info[Y][512 * 3 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 2 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 4 + (415 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 0 + (403 + tid + 256) & 511];
		info1_3 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 4 + (184 + tid + 256) & 511];
		info1_4 = buf_d_channel_info[Y][512 * 0 + (279 + tid + 256) & 511];
		info1_5 = buf_d_channel_info[Y][512 * 3 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid) & 511)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid) & 511)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
		info1_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos1_3];
		info1_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos1_4];
		info1_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos1_5];
	}

	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	symbol0 |= info_symbol0 << 0;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	symbol0 |= info_symbol1 << 1;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	symbol0 |= info_symbol2 << 2;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	symbol0 |= info_symbol3 << 3;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	symbol0 |= info_symbol4 << 4;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	symbol0 |= info_symbol5 << 5;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 0 & 0x01)));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 1 & 0x01)));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 2 & 0x01)));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 3 & 0x01)));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 4 & 0x01)));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 5 & 0x01)));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	symbol1 |= info_symbol0 << 0;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	symbol1 |= info_symbol1 << 1;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	symbol1 |= info_symbol2 << 2;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;

	info_symbol3 = ((unsigned int)(int)(info1_3 - 1)) >> 31;
	symbol1 |= info_symbol3 << 3;
	info1_3 = dLtanh(info1_3);
	info_sum_symbol1 ^= info_symbol3;
	info_sum1 += info1_3;

	info_symbol4 = ((unsigned int)(int)(info1_4 - 1)) >> 31;
	symbol1 |= info_symbol4 << 4;
	info1_4 = dLtanh(info1_4);
	info_sum_symbol1 ^= info_symbol4;
	info_sum1 += info1_4;

	info_symbol5 = ((unsigned int)(int)(info1_5 - 1)) >> 31;
	symbol1 |= info_symbol5 << 5;
	info1_5 = dLtanh(info1_5);
	info_sum_symbol1 ^= info_symbol5;
	info_sum1 += info1_5;

	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 0 & 0x01)));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 1 & 0x01)));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 2 & 0x01)));

	info1_3 = info_sum1 - info1_3;
	info1_3 = dLtanh(info1_3);
	info1_3 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 3 & 0x01)));

	info1_4 = info_sum1 - info1_4;
	info1_4 = dLtanh(info1_4);
	info1_4 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 4 & 0x01)));

	info1_5 = info_sum1 - info1_5;
	info1_5 = dLtanh(info1_5);
	info1_5 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 5 & 0x01)));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;
	d_info_row_2_col[Y].info[3][(row_id + 256) % ROW_LENGTH] = info1_3;
	d_info_row_2_col[Y].info[4][(row_id + 256) % ROW_LENGTH] = info1_4;
	d_info_row_2_col[Y].info[5][(row_id + 256) % ROW_LENGTH] = info1_5;



	/**************************ROW_2********************************/
	row_id = row_offset - I * 4 * 1536 + 2 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0; symbol0 = 0; symbol1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 3 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 0 + (48 + tid) & 511];
		info0_2 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 1 + (7 + tid) & 511];
		info0_3 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 2 + (328 + tid) & 511];
		info0_4 = buf_d_channel_info[Y][512 * 0 + (185 + tid) & 511];
		info0_5 = buf_d_channel_info[Y][512 * 4 + (0 + tid) & 511];

		info0_0 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 3 + (0 + tid + 256) & 511];
		info0_1 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 0 + (48 + tid + 256) & 511];
		info0_2 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 1 + (7 + tid + 256) & 511];
		info0_3 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 2 + (328 + tid + 256) & 511];
		info0_4 = buf_d_channel_info[Y][512 * 0 + (185 + tid + 256) & 511];
		info0_5 = buf_d_channel_info[Y][512 * 4 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid) & 511)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid) & 511)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
		info1_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos1_3];
		info1_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos1_4];
		info1_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos1_5];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	symbol0 |= info_symbol0 << 0;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	symbol0 |= info_symbol1 << 1;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	symbol0 |= info_symbol2 << 2;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	symbol0 |= info_symbol3 << 3;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	symbol0 |= info_symbol4 << 4;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	symbol0 |= info_symbol5 << 5;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 0 & 0x01)));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 1 & 0x01)));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 2 & 0x01)));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 3 & 0x01)));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 4 & 0x01)));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 5 & 0x01)));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	symbol1 |= info_symbol0 << 0;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	symbol1 |= info_symbol1 << 1;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	symbol1 |= info_symbol2 << 2;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;

	info_symbol3 = ((unsigned int)(int)(info1_3 - 1)) >> 31;
	symbol1 |= info_symbol3 << 3;
	info1_3 = dLtanh(info1_3);
	info_sum_symbol1 ^= info_symbol3;
	info_sum1 += info1_3;

	info_symbol4 = ((unsigned int)(int)(info1_4 - 1)) >> 31;
	symbol1 |= info_symbol4 << 4;
	info1_4 = dLtanh(info1_4);
	info_sum_symbol1 ^= info_symbol4;
	info_sum1 += info1_4;

	info_symbol5 = ((unsigned int)(int)(info1_5 - 1)) >> 31;
	symbol1 |= info_symbol5 << 5;
	info1_5 = dLtanh(info1_5);
	info_sum_symbol1 ^= info_symbol5;
	info_sum1 += info1_5;

	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 0 & 0x01)));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 1 & 0x01)));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 2 & 0x01)));

	info1_3 = info_sum1 - info1_3;
	info1_3 = dLtanh(info1_3);
	info1_3 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 3 & 0x01)));

	info1_4 = info_sum1 - info1_4;
	info1_4 = dLtanh(info1_4);
	info1_4 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 4 & 0x01)));

	info1_5 = info_sum1 - info1_5;
	info1_5 = dLtanh(info1_5);
	info1_5 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 5 & 0x01)));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;
	d_info_row_2_col[Y].info[3][(row_id + 256) % ROW_LENGTH] = info1_3;
	d_info_row_2_col[Y].info[4][(row_id + 256) % ROW_LENGTH] = info1_4;
	d_info_row_2_col[Y].info[5][(row_id + 256) % ROW_LENGTH] = info1_5;


	__syncthreads();



	/***************************************************************/
	/***********************COL_UPDATE******************************/
	/***************************************************************/

	/**************************COL_0********************************/
	col_id = col_offset - I * 4 * 2560 + 0 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][0 * 512 + tid];
		info_channel1 = buf_d_channel_info[1 * STREAM_COUNT + Y][0 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}

	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];


	info_sum0 = info0_0 + info0_1;
	info_sum1 = info1_0 + info1_1;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(0 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(0 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
	}



	/**************************COL_1********************************/
	col_id = col_offset - I * 4 * 2560 + 1 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][1 * 512 + tid];
		info_channel1 = buf_d_channel_info[1 * STREAM_COUNT + Y][1 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];


	info_sum0 = info0_0 + info0_1;
	info_sum1 = info1_0 + info1_1;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(1 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(1 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
	}



	/**************************COL_2********************************/
	col_id = col_offset - I * 4 * 2560 + 2 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][2 * 512 + tid];
		info_channel1 = buf_d_channel_info[1 * STREAM_COUNT + Y][2 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];


	info_sum0 = info0_0 + info0_1;
	info_sum1 = info1_0 + info1_1;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(2 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(2 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
	}



	/**************************COL_3********************************/
	col_id = col_offset - I * 4 * 2560 + 3 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][3 * 512 + tid];
		info_channel1 = buf_d_channel_info[1 * STREAM_COUNT + Y][3 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];


	info_sum0 = info0_0 + info0_1;
	info_sum1 = info1_0 + info1_1;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(3 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(3 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
	}



	/**************************COL_4********************************/
	col_id = col_offset - I * 4 * 2560 + 4 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[1 * STREAM_COUNT + Y][4 * 512 + tid];
		info_channel1 = buf_d_channel_info[1 * STREAM_COUNT + Y][4 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];


	info_sum0 = info0_0 + info0_1 + info0_2;
	info_sum1 = info1_0 + info1_1 + info1_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(4 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(4 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
	}


	__syncthreads();
};


__global__ void dUpdate2(int time_count, info_ch *d_channel_info, int *d_decoded_word, INFO_COL *d_info_col_2_row, INFO_ROW *d_info_row_2_col
	, buf_info_ch *buf_d_channel_info, BUF_INFO_COL *buf_d_info_col_2_row, int STREAM_COUNT
	)
{
	int tid = threadIdx.x;
	int I = blockIdx.x;
	int Y = blockIdx.y;

	int row_offset = time_count * 1536;
	int col_offset = (time_count - 3) * 2560;
	int row_id = 0;
	int col_id = 0;
	register int node_pos0_0, node_pos0_1, node_pos0_2, node_pos0_3, node_pos0_4, node_pos0_5;
	register int node_pos1_0, node_pos1_1, node_pos1_2, node_pos1_3, node_pos1_4, node_pos1_5;
	register int node_cpos0, node_cpos1, node_cpos2, node_cpos3, node_cpos4, node_cpos5;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	register float info0_0, info0_1, info0_2, info0_3, info0_4, info0_5;
	register float info1_0, info1_1, info1_2, info1_3, info1_4, info1_5;
	register int info_symbol0, info_symbol1, info_symbol2, info_symbol3, info_symbol4, info_symbol5;
	register float info_sum0, info_sum1;
	register int info_sum_symbol0, info_sum_symbol1;
	register float info_channel0, info_channel1;
	int node_pos_2 = 0;


	/***************************************************************/
	/***********************ROW_UPDATE******************************/
	/***************************************************************/

	/**************************ROW_0********************************/
	row_id = row_offset - I * 4 * 1536 + 0 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 3 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[Y][512 * 1 + (0 + tid) & 511];
		info0_2 = buf_d_channel_info[Y][512 * 4 + (108 + tid) & 511];

		info1_0 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 3 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[Y][512 * 1 + (0 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[Y][512 * 4 + (108 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;

		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
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
	info0_0 *= (1 - 2 * (info_sum_symbol0^info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0^info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0^info_symbol2));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;


	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1^info_symbol0));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1^info_symbol1));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1^info_symbol2));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;



	/**************************ROW_1********************************/
	row_id = row_offset - I * 4 * 1536 + 1 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 4 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[Y][512 * 1 + (0 + tid) & 511];
		info0_2 = buf_d_channel_info[Y][512 * 2 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 4 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[Y][512 * 1 + (0 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[Y][512 * 2 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;

		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
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
	info0_0 *= (1 - 2 * (info_sum_symbol0^info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0^info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0^info_symbol2));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;


	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1^info_symbol0));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1^info_symbol1));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1^info_symbol2));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;



	/**************************ROW_2********************************/
	row_id = row_offset - I * 4 * 1536 + 2 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 0 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[Y][512 * 2 + (0 + tid) & 511];
		info0_2 = buf_d_channel_info[Y][512 * 3 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 0 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[Y][512 * 2 + (0 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[Y][512 * 3 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;

		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
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
	info0_0 *= (1 - 2 * (info_sum_symbol0^info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0^info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0^info_symbol2));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;


	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1^info_symbol0));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1^info_symbol1));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1^info_symbol2));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;


	__syncthreads();




	/***************************************************************/
	/***********************COL_UPDATE******************************/
	/***************************************************************/

	/**************************COL_0********************************/
	col_id = col_offset - I * 4 * 2560 + 0 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][0 * 512 + tid];
		info_channel1 = buf_d_channel_info[2 * STREAM_COUNT + Y][0 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];


	info_sum0 = info0_0 + info0_1 + info0_2;
	info_sum1 = info1_0 + info1_1 + info1_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(0 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(0 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
	}



	/**************************COL_1********************************/
	col_id = col_offset - I * 4 * 2560 + 1 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][1 * 512 + tid];
		info_channel1 = buf_d_channel_info[2 * STREAM_COUNT + Y][1 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];


	info_sum0 = info0_0 + info0_1 + info0_2;
	info_sum1 = info1_0 + info1_1 + info1_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(1 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(1 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
	}


	/**************************COL_2********************************/
	col_id = col_offset - I * 4 * 2560 + 2 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][2 * 512 + tid];
		info_channel1 = buf_d_channel_info[2 * STREAM_COUNT + Y][2 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}

	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];


	info_sum0 = info0_0 + info0_1 + info0_2;
	info_sum1 = info1_0 + info1_1 + info1_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(2 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(2 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
	}



	/**************************COL_3********************************/
	col_id = col_offset - I * 4 * 2560 + 3 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][3 * 512 + tid];
		info_channel1 = buf_d_channel_info[2 * STREAM_COUNT + Y][3 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];


	info_sum0 = info0_0;
	info_sum1 = info1_0;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(3 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(3 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
	}



	/**************************COL_4********************************/
	col_id = col_offset - I * 4 * 2560 + 4 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[2 * STREAM_COUNT + Y][4 * 512 + tid];
		info_channel1 = buf_d_channel_info[2 * STREAM_COUNT + Y][4 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];


	info_sum0 = info0_0;
	info_sum1 = info1_0;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(4 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(4 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
	}


	__syncthreads();
};


__global__ void dUpdate3(int time_count, info_ch *d_channel_info, int *d_decoded_word, INFO_COL *d_info_col_2_row, INFO_ROW *d_info_row_2_col
	, buf_info_ch *buf_d_channel_info, BUF_INFO_COL *buf_d_info_col_2_row, int STREAM_COUNT
	)
{
	int tid = threadIdx.x;
	int I = blockIdx.x;
	int Y = blockIdx.y;

	int row_offset = time_count * 1536;
	int col_offset = (time_count - 3) * 2560;
	int row_id = 0;
	int col_id = 0;
	register int node_pos0_0, node_pos0_1, node_pos0_2, node_pos0_3, node_pos0_4, node_pos0_5;
	register int node_pos1_0, node_pos1_1, node_pos1_2, node_pos1_3, node_pos1_4, node_pos1_5;
	register int node_cpos0, node_cpos1, node_cpos2, node_cpos3, node_cpos4, node_cpos5;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	register float info0_0, info0_1, info0_2, info0_3, info0_4, info0_5;
	register float info1_0, info1_1, info1_2, info1_3, info1_4, info1_5;
	register int info_symbol0, info_symbol1, info_symbol2, info_symbol3, info_symbol4, info_symbol5;
	register float info_sum0, info_sum1;
	register int info_sum_symbol0, info_sum_symbol1;
	register float info_channel0, info_channel1;
	int node_pos_2 = 0;


	/***************************************************************/
	/***********************ROW_UPDATE******************************/
	/***************************************************************/

	/**************************ROW_0********************************/
	row_id = row_offset - I * 4 * 1536 + 0 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 1 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[Y][512 * 3 + (0 + tid) & 511];
		info0_2 = buf_d_channel_info[Y][512 * 4 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 1 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[Y][512 * 3 + (0 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[Y][512 * 4 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;

		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
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
	info0_0 *= (1 - 2 * (info_sum_symbol0^info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0^info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0^info_symbol2));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;


	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1^info_symbol0));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1^info_symbol1));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1^info_symbol2));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;



	/**************************ROW_1********************************/
	row_id = row_offset - I * 4 * 1536 + 1 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 2 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[Y][512 * 1 + (126 + tid) & 511];
		info0_2 = buf_d_channel_info[Y][512 * 2 + (238 + tid) & 511];
		info0_3 = buf_d_channel_info[Y][512 * 3 + (481 + tid) & 511];
		info0_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 0 + (0 + tid) & 511];
		info0_5 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 4 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 2 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[Y][512 * 1 + (126 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[Y][512 * 2 + (238 + tid + 256) & 511];
		info1_3 = buf_d_channel_info[Y][512 * 3 + (481 + tid + 256) & 511];
		info1_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 0 + (0 + tid + 256) & 511];
		info1_5 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 4 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid) & 511)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid) & 511)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
		info1_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos1_3];
		info1_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos1_4];
		info1_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos1_5];
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

	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0^info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0^info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0^info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0^info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0^info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0^info_symbol5));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;

	info_symbol3 = ((unsigned int)(int)(info1_3 - 1)) >> 31;
	info1_3 = dLtanh(info1_3);
	info_sum_symbol1 ^= info_symbol3;
	info_sum1 += info1_3;

	info_symbol4 = ((unsigned int)(int)(info1_4 - 1)) >> 31;
	info1_4 = dLtanh(info1_4);
	info_sum_symbol1 ^= info_symbol4;
	info_sum1 += info1_4;

	info_symbol5 = ((unsigned int)(int)(info1_5 - 1)) >> 31;
	info1_5 = dLtanh(info1_5);
	info_sum_symbol1 ^= info_symbol5;
	info_sum1 += info1_5;

	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1^info_symbol0));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1^info_symbol1));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1^info_symbol2));

	info1_3 = info_sum1 - info1_3;
	info1_3 = dLtanh(info1_3);
	info1_3 *= (1 - 2 * (info_sum_symbol1^info_symbol3));

	info1_4 = info_sum1 - info1_4;
	info1_4 = dLtanh(info1_4);
	info1_4 *= (1 - 2 * (info_sum_symbol1^info_symbol4));

	info1_5 = info_sum1 - info1_5;
	info1_5 = dLtanh(info1_5);
	info1_5 *= (1 - 2 * (info_sum_symbol1^info_symbol5));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;
	d_info_row_2_col[Y].info[3][(row_id + 256) % ROW_LENGTH] = info1_3;
	d_info_row_2_col[Y].info[4][(row_id + 256) % ROW_LENGTH] = info1_4;
	d_info_row_2_col[Y].info[5][(row_id + 256) % ROW_LENGTH] = info1_5;



	/**************************ROW_2********************************/
	row_id = row_offset - I * 4 * 1536 + 2 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 3 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[Y][512 * 2 + (375 + tid) & 511];
		info0_2 = buf_d_channel_info[Y][512 * 3 + (436 + tid) & 511];
		info0_3 = buf_d_channel_info[Y][512 * 4 + (350 + tid) & 511];
		info0_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 1 + (0 + tid) & 511];
		info0_5 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 0 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 3 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[Y][512 * 2 + (375 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[Y][512 * 3 + (436 + tid + 256) & 511];
		info1_3 = buf_d_channel_info[Y][512 * 4 + (350 + tid + 256) & 511];
		info1_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 1 + (0 + tid + 256) & 511];
		info1_5 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 0 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid) & 511)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid) & 511)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
		info1_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos1_3];
		info1_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos1_4];
		info1_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos1_5];
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

	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0^info_symbol0));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0^info_symbol1));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0^info_symbol2));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0^info_symbol3));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0^info_symbol4));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0^info_symbol5));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;

	info_symbol3 = ((unsigned int)(int)(info1_3 - 1)) >> 31;
	info1_3 = dLtanh(info1_3);
	info_sum_symbol1 ^= info_symbol3;
	info_sum1 += info1_3;

	info_symbol4 = ((unsigned int)(int)(info1_4 - 1)) >> 31;
	info1_4 = dLtanh(info1_4);
	info_sum_symbol1 ^= info_symbol4;
	info_sum1 += info1_4;

	info_symbol5 = ((unsigned int)(int)(info1_5 - 1)) >> 31;
	info1_5 = dLtanh(info1_5);
	info_sum_symbol1 ^= info_symbol5;
	info_sum1 += info1_5;

	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1^info_symbol0));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1^info_symbol1));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1^info_symbol2));

	info1_3 = info_sum1 - info1_3;
	info1_3 = dLtanh(info1_3);
	info1_3 *= (1 - 2 * (info_sum_symbol1^info_symbol3));

	info1_4 = info_sum1 - info1_4;
	info1_4 = dLtanh(info1_4);
	info1_4 *= (1 - 2 * (info_sum_symbol1^info_symbol4));

	info1_5 = info_sum1 - info1_5;
	info1_5 = dLtanh(info1_5);
	info1_5 *= (1 - 2 * (info_sum_symbol1^info_symbol5));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;
	d_info_row_2_col[Y].info[3][(row_id + 256) % ROW_LENGTH] = info1_3;
	d_info_row_2_col[Y].info[4][(row_id + 256) % ROW_LENGTH] = info1_4;
	d_info_row_2_col[Y].info[5][(row_id + 256) % ROW_LENGTH] = info1_5;


	__syncthreads();



	/***************************************************************/
	/***********************COL_UPDATE******************************/
	/***************************************************************/

	/**************************COL_0********************************/
	col_id = col_offset - I * 4 * 2560 + 0 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][0 * 512 + tid];
		info_channel1 = buf_d_channel_info[3 * STREAM_COUNT + Y][0 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];


	info_sum0 = info0_0;
	info_sum1 = info1_0;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(0 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(0 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
	}



	/**************************COL_1********************************/
	col_id = col_offset - I * 4 * 2560 + 1 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][1 * 512 + tid];
		info_channel1 = buf_d_channel_info[3 * STREAM_COUNT + Y][1 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];


	info_sum0 = info0_0;
	info_sum1 = info1_0;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(1 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(1 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
	}



	/**************************COL_2********************************/
	col_id = col_offset - I * 4 * 2560 + 2 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][2 * 512 + tid];
		info_channel1 = buf_d_channel_info[3 * STREAM_COUNT + Y][2 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];


	info_sum0 = info0_0 + info0_1 + info0_2;
	info_sum1 = info1_0 + info1_1 + info1_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(2 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(2 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
	}



	/**************************COL_3********************************/
	col_id = col_offset - I * 4 * 2560 + 3 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][3 * 512 + tid];
		info_channel1 = buf_d_channel_info[3 * STREAM_COUNT + Y][3 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];


	info_sum0 = info0_0 + info0_1 + info0_2;
	info_sum1 = info1_0 + info1_1 + info1_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(3 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(3 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
	}



	/**************************COL_4********************************/
	col_id = col_offset - I * 4 * 2560 + 4 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[3 * STREAM_COUNT + Y][4 * 512 + tid];
		info_channel1 = buf_d_channel_info[3 * STREAM_COUNT + Y][4 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];


	info_sum0 = info0_0 + info0_1 + info0_2;
	info_sum1 = info1_0 + info1_1 + info1_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(4 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(4 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
	}


	__syncthreads();
};


__global__ void dUpdate4(int time_count, info_ch *d_channel_info, int *d_decoded_word, INFO_COL *d_info_col_2_row, INFO_ROW *d_info_row_2_col
	, buf_info_ch *buf_d_channel_info, BUF_INFO_COL *buf_d_info_col_2_row, int STREAM_COUNT
	)
{
	int tid = threadIdx.x;
	int I = blockIdx.x;
	int Y = blockIdx.y;

	int row_offset = time_count * 1536;
	int col_offset = (time_count - 3) * 2560;
	int row_id = 0;
	int col_id = 0;
	register int node_pos0_0, node_pos0_1, node_pos0_2, node_pos0_3, node_pos0_4, node_pos0_5;
	register int node_pos1_0, node_pos1_1, node_pos1_2, node_pos1_3, node_pos1_4, node_pos1_5;
	register int node_cpos0, node_cpos1, node_cpos2, node_cpos3, node_cpos4, node_cpos5;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	register float info0_0, info0_1, info0_2, info0_3, info0_4, info0_5;
	register float info1_0, info1_1, info1_2, info1_3, info1_4, info1_5;
	register int info_symbol0, info_symbol1, info_symbol2, info_symbol3, info_symbol4, info_symbol5;
	register int symbol0, symbol1;
	register float info_sum0, info_sum1;
	register int info_sum_symbol0, info_sum_symbol1;
	register float info_channel0, info_channel1;
	int node_pos_2 = 0;


	/***************************************************************/
	/***********************ROW_UPDATE******************************/
	/***************************************************************/

	/**************************ROW_0********************************/
	row_id = row_offset - I * 4 * 1536 + 0 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0; symbol0 = 0; symbol1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[Y][512 * 1 + (263 + tid) & 511];
		info0_1 = buf_d_channel_info[Y][512 * 3 + (219 + tid) & 511];
		info0_2 = buf_d_channel_info[Y][512 * 4 + (16 + tid) & 511];
		info0_3 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 2 + (0 + tid) & 511];
		info0_4 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 1 + (0 + tid) & 511];
		info0_5 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 4 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[Y][512 * 1 + (263 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[Y][512 * 3 + (219 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[Y][512 * 4 + (16 + tid + 256) & 511];
		info1_3 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 2 + (0 + tid + 256) & 511];
		info1_4 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 1 + (0 + tid + 256) & 511];
		info1_5 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 4 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid) & 511)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid) & 511)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
		info1_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos1_3];
		info1_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos1_4];
		info1_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos1_5];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	symbol0 |= info_symbol0 << 0;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	symbol0 |= info_symbol1 << 1;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	symbol0 |= info_symbol2 << 2;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	symbol0 |= info_symbol3 << 3;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	symbol0 |= info_symbol4 << 4;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	symbol0 |= info_symbol5 << 5;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 0 & 0x01)));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 1 & 0x01)));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 2 & 0x01)));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 3 & 0x01)));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 4 & 0x01)));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 5 & 0x01)));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	symbol1 |= info_symbol0 << 0;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	symbol1 |= info_symbol1 << 1;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	symbol1 |= info_symbol2 << 2;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;

	info_symbol3 = ((unsigned int)(int)(info1_3 - 1)) >> 31;
	symbol1 |= info_symbol3 << 3;
	info1_3 = dLtanh(info1_3);
	info_sum_symbol1 ^= info_symbol3;
	info_sum1 += info1_3;

	info_symbol4 = ((unsigned int)(int)(info1_4 - 1)) >> 31;
	symbol1 |= info_symbol4 << 4;
	info1_4 = dLtanh(info1_4);
	info_sum_symbol1 ^= info_symbol4;
	info_sum1 += info1_4;

	info_symbol5 = ((unsigned int)(int)(info1_5 - 1)) >> 31;
	symbol1 |= info_symbol5 << 5;
	info1_5 = dLtanh(info1_5);
	info_sum_symbol1 ^= info_symbol5;
	info_sum1 += info1_5;

	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 0 & 0x01)));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 1 & 0x01)));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 2 & 0x01)));

	info1_3 = info_sum1 - info1_3;
	info1_3 = dLtanh(info1_3);
	info1_3 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 3 & 0x01)));

	info1_4 = info_sum1 - info1_4;
	info1_4 = dLtanh(info1_4);
	info1_4 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 4 & 0x01)));

	info1_5 = info_sum1 - info1_5;
	info1_5 = dLtanh(info1_5);
	info1_5 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 5 & 0x01)));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;
	d_info_row_2_col[Y].info[3][(row_id + 256) % ROW_LENGTH] = info1_3;
	d_info_row_2_col[Y].info[4][(row_id + 256) % ROW_LENGTH] = info1_4;
	d_info_row_2_col[Y].info[5][(row_id + 256) % ROW_LENGTH] = info1_5;



	/**************************ROW_1********************************/
	row_id = row_offset - I * 4 * 1536 + 1 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0; symbol0 = 0; symbol1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[Y][512 * 0 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[Y][512 * 1 + (503 + tid) & 511];
		info0_2 = buf_d_channel_info[Y][512 * 2 + (388 + tid) & 511];
		info0_3 = buf_d_channel_info[Y][512 * 4 + (312 + tid) & 511];
		info0_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 3 + (0 + tid) & 511];
		info0_5 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 2 + (0 + tid) & 511];

		info1_0 = buf_d_channel_info[Y][512 * 0 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[Y][512 * 1 + (503 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[Y][512 * 2 + (388 + tid + 256) & 511];
		info1_3 = buf_d_channel_info[Y][512 * 4 + (312 + tid + 256) & 511];
		info1_4 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 3 + (0 + tid + 256) & 511];
		info1_5 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 2 + (0 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid) & 511)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid) & 511)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
		info1_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos1_3];
		info1_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos1_4];
		info1_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos1_5];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	symbol0 |= info_symbol0 << 0;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	symbol0 |= info_symbol1 << 1;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	symbol0 |= info_symbol2 << 2;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	symbol0 |= info_symbol3 << 3;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	symbol0 |= info_symbol4 << 4;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	symbol0 |= info_symbol5 << 5;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 0 & 0x01)));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 1 & 0x01)));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 2 & 0x01)));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 3 & 0x01)));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 4 & 0x01)));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 5 & 0x01)));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	symbol1 |= info_symbol0 << 0;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	symbol1 |= info_symbol1 << 1;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	symbol1 |= info_symbol2 << 2;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;

	info_symbol3 = ((unsigned int)(int)(info1_3 - 1)) >> 31;
	symbol1 |= info_symbol3 << 3;
	info1_3 = dLtanh(info1_3);
	info_sum_symbol1 ^= info_symbol3;
	info_sum1 += info1_3;

	info_symbol4 = ((unsigned int)(int)(info1_4 - 1)) >> 31;
	symbol1 |= info_symbol4 << 4;
	info1_4 = dLtanh(info1_4);
	info_sum_symbol1 ^= info_symbol4;
	info_sum1 += info1_4;

	info_symbol5 = ((unsigned int)(int)(info1_5 - 1)) >> 31;
	symbol1 |= info_symbol5 << 5;
	info1_5 = dLtanh(info1_5);
	info_sum_symbol1 ^= info_symbol5;
	info_sum1 += info1_5;

	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 0 & 0x01)));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 1 & 0x01)));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 2 & 0x01)));

	info1_3 = info_sum1 - info1_3;
	info1_3 = dLtanh(info1_3);
	info1_3 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 3 & 0x01)));

	info1_4 = info_sum1 - info1_4;
	info1_4 = dLtanh(info1_4);
	info1_4 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 4 & 0x01)));

	info1_5 = info_sum1 - info1_5;
	info1_5 = dLtanh(info1_5);
	info1_5 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 5 & 0x01)));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;
	d_info_row_2_col[Y].info[3][(row_id + 256) % ROW_LENGTH] = info1_3;
	d_info_row_2_col[Y].info[4][(row_id + 256) % ROW_LENGTH] = info1_4;
	d_info_row_2_col[Y].info[5][(row_id + 256) % ROW_LENGTH] = info1_5;



	/**************************ROW_2********************************/
	row_id = row_offset - I * 4 * 1536 + 2 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0; symbol0 = 0; symbol1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;

	block_id = (row_id % 6144) >> 9;
	offset = ((row_id / 1536) - ROW_OFFSET) * 2560;

	if (I == 0)
	{
		info0_0 = buf_d_channel_info[Y][512 * 1 + (0 + tid) & 511];
		info0_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 0 + (0 + tid) & 511];
		info0_2 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 1 + (96 + tid) & 511];
		info0_3 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 2 + (28 + tid) & 511];
		info0_4 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 2 + (59 + tid) & 511];
		info0_5 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 3 + (225 + tid) & 511];

		info1_0 = buf_d_channel_info[Y][512 * 1 + (0 + tid + 256) & 511];
		info1_1 = buf_d_channel_info[1 * STREAM_COUNT + Y][512 * 0 + (0 + tid + 256) & 511];
		info1_2 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 1 + (96 + tid + 256) & 511];
		info1_3 = buf_d_channel_info[2 * STREAM_COUNT + Y][512 * 2 + (28 + tid + 256) & 511];
		info1_4 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 2 + (59 + tid + 256) & 511];
		info1_5 = buf_d_channel_info[3 * STREAM_COUNT + Y][512 * 3 + (225 + tid + 256) & 511];
	}
	else
	{
		node_pos0_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid) & 511)) % COL_LENGTH;
		node_pos0_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid) & 511)) % COL_LENGTH;
		node_pos0_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid) & 511)) % COL_LENGTH;
		node_pos0_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid) & 511)) % COL_LENGTH;
		node_pos0_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid) & 511)) % COL_LENGTH;
		node_pos0_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid) & 511)) % COL_LENGTH;

		node_pos1_0 = (offset + d_matrix_node_c[block_id].order[0].col * 512 + ((d_matrix_node_c[block_id].order[0].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_1 = (offset + d_matrix_node_c[block_id].order[1].col * 512 + ((d_matrix_node_c[block_id].order[1].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_2 = (offset + d_matrix_node_c[block_id].order[2].col * 512 + ((d_matrix_node_c[block_id].order[2].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_3 = (offset + d_matrix_node_c[block_id].order[3].col * 512 + ((d_matrix_node_c[block_id].order[3].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_4 = (offset + d_matrix_node_c[block_id].order[4].col * 512 + ((d_matrix_node_c[block_id].order[4].number + tid + 256) & 511)) % COL_LENGTH;
		node_pos1_5 = (offset + d_matrix_node_c[block_id].order[5].col * 512 + ((d_matrix_node_c[block_id].order[5].number + tid + 256) & 511)) % COL_LENGTH;

		node_cpos0 = (d_matrix_node_c[block_id].order[0].col_2_row) % COL_LENGTH;
		node_cpos1 = (d_matrix_node_c[block_id].order[1].col_2_row) % COL_LENGTH;
		node_cpos2 = (d_matrix_node_c[block_id].order[2].col_2_row) % COL_LENGTH;
		node_cpos3 = (d_matrix_node_c[block_id].order[3].col_2_row) % COL_LENGTH;
		node_cpos4 = (d_matrix_node_c[block_id].order[4].col_2_row) % COL_LENGTH;
		node_cpos5 = (d_matrix_node_c[block_id].order[5].col_2_row) % COL_LENGTH;


		info0_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos0_0];
		info0_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos0_1];
		info0_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos0_2];
		info0_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos0_3];
		info0_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos0_4];
		info0_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos0_5];

		info1_0 = d_info_col_2_row[Y].info[node_cpos0][node_pos1_0];
		info1_1 = d_info_col_2_row[Y].info[node_cpos1][node_pos1_1];
		info1_2 = d_info_col_2_row[Y].info[node_cpos2][node_pos1_2];
		info1_3 = d_info_col_2_row[Y].info[node_cpos3][node_pos1_3];
		info1_4 = d_info_col_2_row[Y].info[node_cpos4][node_pos1_4];
		info1_5 = d_info_col_2_row[Y].info[node_cpos5][node_pos1_5];
	}


	info_symbol0 = ((unsigned int)(int)(info0_0 - 1)) >> 31;
	symbol0 |= info_symbol0 << 0;
	info0_0 = dLtanh(info0_0);
	info_sum_symbol0 ^= info_symbol0;
	info_sum0 += info0_0;

	info_symbol1 = ((unsigned int)(int)(info0_1 - 1)) >> 31;
	symbol0 |= info_symbol1 << 1;
	info0_1 = dLtanh(info0_1);
	info_sum_symbol0 ^= info_symbol1;
	info_sum0 += info0_1;

	info_symbol2 = ((unsigned int)(int)(info0_2 - 1)) >> 31;
	symbol0 |= info_symbol2 << 2;
	info0_2 = dLtanh(info0_2);
	info_sum_symbol0 ^= info_symbol2;
	info_sum0 += info0_2;

	info_symbol3 = ((unsigned int)(int)(info0_3 - 1)) >> 31;
	symbol0 |= info_symbol3 << 3;
	info0_3 = dLtanh(info0_3);
	info_sum_symbol0 ^= info_symbol3;
	info_sum0 += info0_3;

	info_symbol4 = ((unsigned int)(int)(info0_4 - 1)) >> 31;
	symbol0 |= info_symbol4 << 4;
	info0_4 = dLtanh(info0_4);
	info_sum_symbol0 ^= info_symbol4;
	info_sum0 += info0_4;

	info_symbol5 = ((unsigned int)(int)(info0_5 - 1)) >> 31;
	symbol0 |= info_symbol5 << 5;
	info0_5 = dLtanh(info0_5);
	info_sum_symbol0 ^= info_symbol5;
	info_sum0 += info0_5;

	info0_0 = info_sum0 - info0_0;
	info0_0 = dLtanh(info0_0);
	info0_0 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 0 & 0x01)));

	info0_1 = info_sum0 - info0_1;
	info0_1 = dLtanh(info0_1);
	info0_1 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 1 & 0x01)));

	info0_2 = info_sum0 - info0_2;
	info0_2 = dLtanh(info0_2);
	info0_2 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 2 & 0x01)));

	info0_3 = info_sum0 - info0_3;
	info0_3 = dLtanh(info0_3);
	info0_3 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 3 & 0x01)));

	info0_4 = info_sum0 - info0_4;
	info0_4 = dLtanh(info0_4);
	info0_4 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 4 & 0x01)));

	info0_5 = info_sum0 - info0_5;
	info0_5 = dLtanh(info0_5);
	info0_5 *= (1 - 2 * (info_sum_symbol0 ^ (symbol0 >> 5 & 0x01)));


	info_symbol0 = ((unsigned int)(int)(info1_0 - 1)) >> 31;
	symbol1 |= info_symbol0 << 0;
	info1_0 = dLtanh(info1_0);
	info_sum_symbol1 ^= info_symbol0;
	info_sum1 += info1_0;

	info_symbol1 = ((unsigned int)(int)(info1_1 - 1)) >> 31;
	symbol1 |= info_symbol1 << 1;
	info1_1 = dLtanh(info1_1);
	info_sum_symbol1 ^= info_symbol1;
	info_sum1 += info1_1;

	info_symbol2 = ((unsigned int)(int)(info1_2 - 1)) >> 31;
	symbol1 |= info_symbol2 << 2;
	info1_2 = dLtanh(info1_2);
	info_sum_symbol1 ^= info_symbol2;
	info_sum1 += info1_2;

	info_symbol3 = ((unsigned int)(int)(info1_3 - 1)) >> 31;
	symbol1 |= info_symbol3 << 3;
	info1_3 = dLtanh(info1_3);
	info_sum_symbol1 ^= info_symbol3;
	info_sum1 += info1_3;

	info_symbol4 = ((unsigned int)(int)(info1_4 - 1)) >> 31;
	symbol1 |= info_symbol4 << 4;
	info1_4 = dLtanh(info1_4);
	info_sum_symbol1 ^= info_symbol4;
	info_sum1 += info1_4;

	info_symbol5 = ((unsigned int)(int)(info1_5 - 1)) >> 31;
	symbol1 |= info_symbol5 << 5;
	info1_5 = dLtanh(info1_5);
	info_sum_symbol1 ^= info_symbol5;
	info_sum1 += info1_5;

	info1_0 = info_sum1 - info1_0;
	info1_0 = dLtanh(info1_0);
	info1_0 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 0 & 0x01)));

	info1_1 = info_sum1 - info1_1;
	info1_1 = dLtanh(info1_1);
	info1_1 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 1 & 0x01)));

	info1_2 = info_sum1 - info1_2;
	info1_2 = dLtanh(info1_2);
	info1_2 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 2 & 0x01)));

	info1_3 = info_sum1 - info1_3;
	info1_3 = dLtanh(info1_3);
	info1_3 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 3 & 0x01)));

	info1_4 = info_sum1 - info1_4;
	info1_4 = dLtanh(info1_4);
	info1_4 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 4 & 0x01)));

	info1_5 = info_sum1 - info1_5;
	info1_5 = dLtanh(info1_5);
	info1_5 *= (1 - 2 * (info_sum_symbol1 ^ (symbol1 >> 5 & 0x01)));


	d_info_row_2_col[Y].info[0][row_id%ROW_LENGTH] = info0_0;
	d_info_row_2_col[Y].info[1][row_id%ROW_LENGTH] = info0_1;
	d_info_row_2_col[Y].info[2][row_id%ROW_LENGTH] = info0_2;
	d_info_row_2_col[Y].info[3][row_id%ROW_LENGTH] = info0_3;
	d_info_row_2_col[Y].info[4][row_id%ROW_LENGTH] = info0_4;
	d_info_row_2_col[Y].info[5][row_id%ROW_LENGTH] = info0_5;


	d_info_row_2_col[Y].info[0][(row_id + 256) % ROW_LENGTH] = info1_0;
	d_info_row_2_col[Y].info[1][(row_id + 256) % ROW_LENGTH] = info1_1;
	d_info_row_2_col[Y].info[2][(row_id + 256) % ROW_LENGTH] = info1_2;
	d_info_row_2_col[Y].info[3][(row_id + 256) % ROW_LENGTH] = info1_3;
	d_info_row_2_col[Y].info[4][(row_id + 256) % ROW_LENGTH] = info1_4;
	d_info_row_2_col[Y].info[5][(row_id + 256) % ROW_LENGTH] = info1_5;


	__syncthreads();



	/***************************************************************/
	/***********************COL_UPDATE******************************/
	/***************************************************************/

	/**************************COL_0********************************/
	col_id = col_offset - I * 4 * 2560 + 0 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[Y][0 * 512 + tid];
		info_channel1 = buf_d_channel_info[Y][0 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];


	info_sum0 = info0_0 + info0_1 + info0_2;
	info_sum1 = info1_0 + info1_1 + info1_2;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(0 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(0 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
	}



	/**************************COL_1********************************/
	col_id = col_offset - I * 4 * 2560 + 1 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[Y][1 * 512 + tid];
		info_channel1 = buf_d_channel_info[Y][1 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 512 + ((tid - d_matrix_node_c[block_id].order[3].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_4 = (offset + d_matrix_node_v[block_id].order[4].row * 512 + ((tid - d_matrix_node_c[block_id].order[4].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_5 = (offset + d_matrix_node_v[block_id].order[5].row * 512 + ((tid - d_matrix_node_c[block_id].order[5].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_3 = (offset + d_matrix_node_v[block_id].order[3].row * 512 + ((tid - d_matrix_node_c[block_id].order[3].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_4 = (offset + d_matrix_node_v[block_id].order[4].row * 512 + ((tid - d_matrix_node_c[block_id].order[4].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_5 = (offset + d_matrix_node_v[block_id].order[5].row * 512 + ((tid - d_matrix_node_c[block_id].order[5].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;
	node_cpos4 = (d_matrix_node_v[block_id].order[4].row_2_col) % ROW_LENGTH;
	node_cpos5 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];
	info0_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos0_4];
	info0_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos0_5];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];
	info1_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos1_3];
	info1_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos1_4];
	info1_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos1_5];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3 + info0_4 + info0_5;
	info_sum1 = info1_0 + info1_1 + info1_2 + info1_3 + info1_4 + info1_5;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(1 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(1 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;
		info0_4 = info_sum0 - info0_4 + info_channel0;
		info0_5 = info_sum0 - info0_5 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;
		info1_3 = info_sum1 - info1_3 + info_channel1;
		info1_4 = info_sum1 - info1_4 + info_channel1;
		info1_5 = info_sum1 - info1_5 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
		d_info_col_2_row[Y].info[4][col_id%COL_LENGTH] = info0_4;
		d_info_col_2_row[Y].info[5][col_id%COL_LENGTH] = info0_5;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
		d_info_col_2_row[Y].info[3][(col_id + 256) % COL_LENGTH] = info1_3;
		d_info_col_2_row[Y].info[4][(col_id + 256) % COL_LENGTH] = info1_4;
		d_info_col_2_row[Y].info[5][(col_id + 256) % COL_LENGTH] = info1_5;
	}



	/**************************COL_2********************************/
	col_id = col_offset - I * 4 * 2560 + 2 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[Y][2 * 512 + tid];
		info_channel1 = buf_d_channel_info[Y][2 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 512 + ((tid - d_matrix_node_c[block_id].order[3].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_4 = (offset + d_matrix_node_v[block_id].order[4].row * 512 + ((tid - d_matrix_node_c[block_id].order[4].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_5 = (offset + d_matrix_node_v[block_id].order[5].row * 512 + ((tid - d_matrix_node_c[block_id].order[5].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_3 = (offset + d_matrix_node_v[block_id].order[3].row * 512 + ((tid - d_matrix_node_c[block_id].order[3].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_4 = (offset + d_matrix_node_v[block_id].order[4].row * 512 + ((tid - d_matrix_node_c[block_id].order[4].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_5 = (offset + d_matrix_node_v[block_id].order[5].row * 512 + ((tid - d_matrix_node_c[block_id].order[5].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;
	node_cpos4 = (d_matrix_node_v[block_id].order[4].row_2_col) % ROW_LENGTH;
	node_cpos5 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];
	info0_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos0_4];
	info0_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos0_5];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];
	info1_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos1_3];
	info1_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos1_4];
	info1_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos1_5];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3 + info0_4 + info0_5;
	info_sum1 = info1_0 + info1_1 + info1_2 + info1_3 + info1_4 + info1_5;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(2 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(2 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;
		info0_4 = info_sum0 - info0_4 + info_channel0;
		info0_5 = info_sum0 - info0_5 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;
		info1_3 = info_sum1 - info1_3 + info_channel1;
		info1_4 = info_sum1 - info1_4 + info_channel1;
		info1_5 = info_sum1 - info1_5 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
		d_info_col_2_row[Y].info[4][col_id%COL_LENGTH] = info0_4;
		d_info_col_2_row[Y].info[5][col_id%COL_LENGTH] = info0_5;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
		d_info_col_2_row[Y].info[3][(col_id + 256) % COL_LENGTH] = info1_3;
		d_info_col_2_row[Y].info[4][(col_id + 256) % COL_LENGTH] = info1_4;
		d_info_col_2_row[Y].info[5][(col_id + 256) % COL_LENGTH] = info1_5;
	}



	/**************************COL_3********************************/
	col_id = col_offset - I * 4 * 2560 + 3 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[Y][3 * 512 + tid];
		info_channel1 = buf_d_channel_info[Y][3 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 512 + ((tid - d_matrix_node_c[block_id].order[3].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_4 = (offset + d_matrix_node_v[block_id].order[4].row * 512 + ((tid - d_matrix_node_c[block_id].order[4].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_5 = (offset + d_matrix_node_v[block_id].order[5].row * 512 + ((tid - d_matrix_node_c[block_id].order[5].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_3 = (offset + d_matrix_node_v[block_id].order[3].row * 512 + ((tid - d_matrix_node_c[block_id].order[3].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_4 = (offset + d_matrix_node_v[block_id].order[4].row * 512 + ((tid - d_matrix_node_c[block_id].order[4].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_5 = (offset + d_matrix_node_v[block_id].order[5].row * 512 + ((tid - d_matrix_node_c[block_id].order[5].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;
	node_cpos4 = (d_matrix_node_v[block_id].order[4].row_2_col) % ROW_LENGTH;
	node_cpos5 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];
	info0_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos0_4];
	info0_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos0_5];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];
	info1_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos1_3];
	info1_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos1_4];
	info1_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos1_5];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3 + info0_4 + info0_5;
	info_sum1 = info1_0 + info1_1 + info1_2 + info1_3 + info1_4 + info1_5;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(3 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(3 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;
		info0_4 = info_sum0 - info0_4 + info_channel0;
		info0_5 = info_sum0 - info0_5 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;
		info1_3 = info_sum1 - info1_3 + info_channel1;
		info1_4 = info_sum1 - info1_4 + info_channel1;
		info1_5 = info_sum1 - info1_5 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
		d_info_col_2_row[Y].info[4][col_id%COL_LENGTH] = info0_4;
		d_info_col_2_row[Y].info[5][col_id%COL_LENGTH] = info0_5;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
		d_info_col_2_row[Y].info[3][(col_id + 256) % COL_LENGTH] = info1_3;
		d_info_col_2_row[Y].info[4][(col_id + 256) % COL_LENGTH] = info1_4;
		d_info_col_2_row[Y].info[5][(col_id + 256) % COL_LENGTH] = info1_5;
	}



	/**************************COL_4********************************/
	col_id = col_offset - I * 4 * 2560 + 4 * 512 + tid;
	info_sum0 = 0; info_sum1 = 0;
	info_sum_symbol0 = 0; info_sum_symbol1 = 0;
	info_channel0 = 0; info_channel1 = 0;

	block_id = (col_id % 10240) >> 9;
	offset = (col_id / 2560) * 1536;

	if (I == 0)
	{
		info_channel0 = buf_d_channel_info[Y][4 * 512 + tid];
		info_channel1 = buf_d_channel_info[Y][4 * 512 + tid + 256];

		d_channel_info[Y][col_id%COL_LENGTH] = info_channel0;
		d_channel_info[Y][(col_id + 256) % COL_LENGTH] = info_channel1;
	}
	else
	{
		info_channel0 = d_channel_info[Y][col_id%COL_LENGTH];
		info_channel1 = d_channel_info[Y][(col_id + 256) % COL_LENGTH];
	}


	node_pos0_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_3 = (offset + d_matrix_node_v[block_id].order[3].row * 512 + ((tid - d_matrix_node_c[block_id].order[3].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_4 = (offset + d_matrix_node_v[block_id].order[4].row * 512 + ((tid - d_matrix_node_c[block_id].order[4].number + 512) & 511)) % ROW_LENGTH;
	node_pos0_5 = (offset + d_matrix_node_v[block_id].order[5].row * 512 + ((tid - d_matrix_node_c[block_id].order[5].number + 512) & 511)) % ROW_LENGTH;

	node_pos1_0 = (offset + d_matrix_node_v[block_id].order[0].row * 512 + ((tid - d_matrix_node_c[block_id].order[0].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_1 = (offset + d_matrix_node_v[block_id].order[1].row * 512 + ((tid - d_matrix_node_c[block_id].order[1].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_2 = (offset + d_matrix_node_v[block_id].order[2].row * 512 + ((tid - d_matrix_node_c[block_id].order[2].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_3 = (offset + d_matrix_node_v[block_id].order[3].row * 512 + ((tid - d_matrix_node_c[block_id].order[3].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_4 = (offset + d_matrix_node_v[block_id].order[4].row * 512 + ((tid - d_matrix_node_c[block_id].order[4].number + 512 + 256) & 511)) % ROW_LENGTH;
	node_pos1_5 = (offset + d_matrix_node_v[block_id].order[5].row * 512 + ((tid - d_matrix_node_c[block_id].order[5].number + 512 + 256) & 511)) % ROW_LENGTH;

	node_cpos0 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;
	node_cpos1 = (d_matrix_node_v[block_id].order[1].row_2_col) % ROW_LENGTH;
	node_cpos2 = (d_matrix_node_v[block_id].order[2].row_2_col) % ROW_LENGTH;
	node_cpos3 = (d_matrix_node_v[block_id].order[3].row_2_col) % ROW_LENGTH;
	node_cpos4 = (d_matrix_node_v[block_id].order[4].row_2_col) % ROW_LENGTH;
	node_cpos5 = (d_matrix_node_v[block_id].order[0].row_2_col) % ROW_LENGTH;


	info0_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos0_0];
	info0_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos0_1];
	info0_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos0_2];
	info0_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos0_3];
	info0_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos0_4];
	info0_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos0_5];

	info1_0 = d_info_row_2_col[Y].info[node_cpos0][node_pos1_0];
	info1_1 = d_info_row_2_col[Y].info[node_cpos1][node_pos1_1];
	info1_2 = d_info_row_2_col[Y].info[node_cpos2][node_pos1_2];
	info1_3 = d_info_row_2_col[Y].info[node_cpos3][node_pos1_3];
	info1_4 = d_info_row_2_col[Y].info[node_cpos4][node_pos1_4];
	info1_5 = d_info_row_2_col[Y].info[node_cpos5][node_pos1_5];


	info_sum0 = info0_0 + info0_1 + info0_2 + info0_3 + info0_4 + info0_5;
	info_sum1 = info1_0 + info1_1 + info1_2 + info1_3 + info1_4 + info1_5;


	if (blockIdx.x == (ITERATE_TIME - 1))
	{
		int decision = 0;
		info_sum0 += info_channel0;
		decision = ((unsigned int)(int)(info_sum0 - 1)) >> 31;
		d_decoded_word[(4 * 512 + tid) % 2560 + Y * 2560] = decision;

		decision = 0;
		info_sum1 += info_channel1;
		decision = ((unsigned int)(int)(info_sum1 - 1)) >> 31;
		d_decoded_word[(4 * 512 + tid + 256) % 2560 + Y * 2560] = decision;
	}
	else
	{
		info0_0 = info_sum0 - info0_0 + info_channel0;
		info0_1 = info_sum0 - info0_1 + info_channel0;
		info0_2 = info_sum0 - info0_2 + info_channel0;
		info0_3 = info_sum0 - info0_3 + info_channel0;
		info0_4 = info_sum0 - info0_4 + info_channel0;
		info0_5 = info_sum0 - info0_5 + info_channel0;

		info1_0 = info_sum1 - info1_0 + info_channel1;
		info1_1 = info_sum1 - info1_1 + info_channel1;
		info1_2 = info_sum1 - info1_2 + info_channel1;
		info1_3 = info_sum1 - info1_3 + info_channel1;
		info1_4 = info_sum1 - info1_4 + info_channel1;
		info1_5 = info_sum1 - info1_5 + info_channel1;

		d_info_col_2_row[Y].info[0][col_id%COL_LENGTH] = info0_0;
		d_info_col_2_row[Y].info[1][col_id%COL_LENGTH] = info0_1;
		d_info_col_2_row[Y].info[2][col_id%COL_LENGTH] = info0_2;
		d_info_col_2_row[Y].info[3][col_id%COL_LENGTH] = info0_3;
		d_info_col_2_row[Y].info[4][col_id%COL_LENGTH] = info0_4;
		d_info_col_2_row[Y].info[5][col_id%COL_LENGTH] = info0_5;

		d_info_col_2_row[Y].info[0][(col_id + 256) % COL_LENGTH] = info1_0;
		d_info_col_2_row[Y].info[1][(col_id + 256) % COL_LENGTH] = info1_1;
		d_info_col_2_row[Y].info[2][(col_id + 256) % COL_LENGTH] = info1_2;
		d_info_col_2_row[Y].info[3][(col_id + 256) % COL_LENGTH] = info1_3;
		d_info_col_2_row[Y].info[4][(col_id + 256) % COL_LENGTH] = info1_4;
		d_info_col_2_row[Y].info[5][(col_id + 256) % COL_LENGTH] = info1_5;
	}


	__syncthreads();
};


#endif
