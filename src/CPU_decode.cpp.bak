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

#include "CPU_decode.h"

extern matrix_check_node h_matrix_node_c[BLOCK_NUM_ROW];
extern matrix_variable_node h_matrix_node_v[BLOCK_NUM_COL];


extern "C"
int countError(int *word)
{
	int error_count = 0;
	for (int i = 0; i<BLOCK_SIZE; i++)
	{
		if (word[i] != 0)
		{
			error_count++;
		}
	}
	return error_count;
}

extern "C"
float dev_ltanh(float info)
{
	float y, x_abs, x_e;

	x_abs = (float)fabs(info);

	if (x_abs < MINLOG)
	{
		y = MAXLOG;
	}
	else if (x_abs > MAXLOG)
	{
		y = 0.0;
	}
	else
	{
		x_e = (float)exp(-x_abs);
		y = (float)log((1 + x_e) / (1 - x_e));

		if (y > MAXLOG)
		{
			y = MAXLOG;
		}
		else if (y < -MAXLOG)
		{
			y = MAXLOG;
		}
	}
	return (y);
}

#ifdef CODE1
extern "C"
void fun_matrix()
{
	// checkMatrix
	int row_deg[12] = {3,3,3, 3,6,6,
		6,6,6, 6,6,6}; //一共12个行块

	check_node c_node [72] =
	{	{3,0,0},	{11,0,0},	{14,1,108},	{3,0,0},	{11,0,0},	{14,1,108},
	{4,0,0},	{11,1,0},	{12,1,0},   {4,0,0},	{11,1,0},	{12,1,0},
	{5,0,0},	{12,2,0},	{13,1,0},   {5,0,0},	{12,2,0},	{13,1,0},
	{1,0,0},	{8,2,0},	{9,2,0},    {1,0,0},	{8,2,0},	{9,2,0},
	{2,2,0},	{6,2,126},	{7,3,238},	{8,3,481},	{10,0,0},	{14,0,0},
	{3,2,0},	{7,4,375},	{8,4,436},	{9,3,350},	{11,0,0},	{15,0,0},
	{1,3,263},	{3,5,219},	{4,4,16},	{7,0,0},	{11,0,0},	{19,0,0},
	{0,2,0},	{1,4,503},	{2,5,388},	{4,5,312},	{8,0,0},	{12,0,0},
	{1,5,0},	{5,1,0},	{11,1,96},	{12,1,28},	{17,0,59},	{18,0,225},
	{1,1,0},	{4,1,84},	{7,2,260},	{13,1,318},	{14,1,382},	{17,0,0},
	{2,1,0},	{4,2,415},	{5,1,403},	{14,2,184},	{15,0,279},	{18,0,0},
	{3,1,0},	{5,2,48},	{6,2,7},	{12,1,328},	{15,1,185},	{19,0,0}
	};  //不够最大行重进行填充


	for (int i=0;i<12;i++)
	{
		h_matrix_node_c[i].row_deg = row_deg[i];
		for (int j=0;j<6;j++)
		{
			h_matrix_node_c[i].order[j] = c_node[i*6 + j];
		}
	}

	//variableMatrix
	int col_deg[20] = {2,2,2,2,3, 3,3,3,1,1,
		1,1,3,3,3, 3,6,6,6,6}; //一共20个列块

	variable_node v_node [120] =
	{
		{4,4,0},  {8,1,0},  {4,4,0},  {8,1,0},  {4,4,0},  {8,1,0},
		{5,4,0},  {9,0,0},  {5,4,0},  {9,0,0},  {5,4,0},  {9,0,0},
		{6,3,0},  {10,0,0},  {6,3,0},  {10,0,0},  {6,3,0},  {10,0,0},
		{7,4,0},  {11,0,0},  {7,3,0},  {11,0,0},  {7,3,0},  {11,0,0},
		{4,5,0},  {9,1,84},  {10,1,415},  {4,5,0},  {9,1,84},  {10,1,415},
		{2,5,0},  {7,2,403},  {8,1,48},  {2,5,0},  {7,2,403},  {8,1,48},
		{3,4,0},  {5,2,96},  {8,2,7},  {3,4,0},  {5,2,96},  {8,2,7},
		{4,5,0},  {5,3,28},  {6,2,260},  {4,4,0},  {5,3,28},  {6,2,260},
		{9,0,0},  {9,0,0},  {9,0,0},  {9,0,0},  {9,0,0},  {9,0,0},
		{10,0,0},  {10,0,0},  {10,0,0},  {10,0,0},  {10,0,0},  {10,0,0},
		{8,0,0},  {8,0,0},  {8,0,0},  {8,0,0},  {8,0,0},  {8,0,0},
		{9,0,0},  {9,0,0},  {9,0,0},  {9,0,0},  {9,0,0},  {9,0,0},
		{2,4,59},  {5,3,328},  {10,0,0},  {2,4,59},  {5,3,328},  {10,0,0},
		{2,5,225},  {3,3,318},  {11,0,0},  {2,5,225},  {3,3,318},  {11,0,0},
		{0,5,0},  {3,4,382},  {4,3,184},  {0,5,0},  {3,4,382},  {4,3,184},
		{1,4,279},  {2,4,185},  {10,0,0},  {2,4,185},  {1,4,279},  {2,4,185},
		{3,1,0},  {4,1,0},  {7,1,126},  {9,0,263},  {10,1,503},  {11,0,0},
		{0,5,0},  {4,2,0},  {5,1,0},  {7,2,238},  {8,1,375},  {10,2,388},
		{1,5,0},  {5,2,0},  {6,1,0},  {7,3,481},  {8,2,436},  {9,1,219},
		{2,5,0},  {3,2,108},  {6,2,0},  {8,3,350},  {9,2,16},  {10,3,312}
	}; //不够最大列重进行填充


	for (int i=0;i<20;i++)
	{
		h_matrix_node_v[i].col_deg = col_deg[i];
		for (int j=0;j<6;j++)
		{
			h_matrix_node_v[i].order[j] = v_node[i*6 + j];
		}
	}
}

extern "C"
void update(int time_count,float *h_channel_info,int **h_decoded_word,INFO_COL *h_info_col_2_row,INFO_ROW *h_info_row_2_col)
{
	int row_offset = time_count*1536; //offset of the first row
	int col_offset = (time_count-3)*2560; //offset of the first col
	int row_id = 0;
	int col_id = 0;
	node_position node[6];
	int relative_id = 0;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	float info[6] = {0,0,0,0,0,0}; //info = channel col to row info
	int info_symbol[6] = {0,0,0,0,0,0}; //positive = 0, negative = 1
	float info_sum = 0;
	int info_sum_symbol = 0;
	float info_channel = 0;

	//Check Update
	for (int I=0;I<ITERATE_TIME;I++) //20 processor
	{
		for (int r=0;r<CHECK_SIZE;r++) //r Row
		{
			row_id = row_offset - I*4*1536 + r;
			i = 0;
			info_sum = 0;
			info_sum_symbol = 0;

			if (row_id >=0)
			{
				relative_id = row_id%512;
				block_id = (row_id%6144)/512;
				deg = h_matrix_node_c[block_id].row_deg;
				offset = ((row_id/1536) - ROW_OFFSET)*2560;
				while(i<deg)
				{
					node[i].pos = h_matrix_node_c[block_id].order[i].col;
					node[i].cpos = h_matrix_node_c[block_id].order[i].col_2_row;
					number = h_matrix_node_c[block_id].order[i].number;
					node[i].pos = offset + node[i].pos*512 + (number + relative_id)%512;
					node[i].pos = node[i].pos%COL_LENGTH; //cycle to store
					//	printf("pos %d  cpos %d \n",node[i].pos,node[i].cpos);

					//info[i] = channel_col_2_row_info[node[i].pos][node[i].cpos];

					if (node[i].pos >= 0)
					{
						info[i] = h_info_col_2_row->info[node[i].cpos][node[i].pos];
					}
					else
					{
						info[i] = 0;
					}

					i++;
				}

				i=0;
				while(i<deg)
				{
					unsigned int* pa = (unsigned int *)(&info[i]);
					info_symbol[i] = *pa>>31;
					if (info[i] != 0)
					{
						info[i] = dev_ltanh(info[i]);
						info_sum += info[i];
						info_sum_symbol ^= info_symbol[i];
					}
					else
					{
						info[i] = 0;
					}
					i++;
				}

				i=0;
				while (i<deg)
				{
					if (info[i] != 0)
					{
						info[i] = info_sum - info[i];// minus itself
						info[i] = dev_ltanh(info[i]);
						unsigned int* pa = (unsigned int *)(&info[i]);
						*pa ^= ((info_sum_symbol^info_symbol[i])<<31);

						//store the info into row_2_col
						//channel_row_2_col_info[id%ROW_LENGTH][i] = info[i];
						h_info_row_2_col->info[i][row_id%ROW_LENGTH] = info[i];
					}
					i++;
				}

			}
		}
	}




	//Variable Update


	for (int I=0;I<ITERATE_TIME;I++)
	{
		for (int c=0;c<BLOCK_SIZE;c++)
		{
			col_id = col_offset - I*4*2560 + c;
			i = 0;
			info_sum = 0;
			info_sum_symbol = 0;
			info_channel = 0;

			if (col_id >= 0)
			{
				relative_id = col_id%512;
				block_id = (col_id%10240)/512;
				deg = h_matrix_node_v[block_id].col_deg;
				number = 0;
				offset = (col_id/2560)*1536;

				info_channel = h_channel_info[col_id%COL_LENGTH];

				while(i<deg)
				{
					node[i].pos = h_matrix_node_v[block_id].order[i].row;
					node[i].cpos = h_matrix_node_v[block_id].order[i].row_2_col;
					number = h_matrix_node_v[block_id].order[i].number;
					node[i].pos = offset + node[i].pos*512 + (relative_id - number + 512)%512;
					node[i].pos = node[i].pos%ROW_LENGTH; //cycle to store
					//		printf("pos %d  cpos %d \n",node[i].pos,node[i].cpos);

					//info[i] = channel_row_2_col_info[node[i].pos][node[i].cpos];
					info[i] = h_info_row_2_col->info[node[i].cpos][node[i].pos];

					i++;
				}

				i=0;
				while(i<deg)
				{
					info_sum += info[i];
					i++;
				}

				i=0;
				while(i<deg)
				{
					float test_info = info[i];
					float test_info_sum = info_sum;
					info[i] = info_sum - info[i] + info_channel;

					h_info_col_2_row->info[i][col_id%COL_LENGTH]= info[i];
					i++;
				}


				//make decision, I = 19
				if ( I == (ITERATE_TIME - 1))
				{
					int decision = 0;

					info_sum += info_channel;
					if (info_sum<0)
					{
						decision = 1;
					}
					else
					{
						decision = 0;
					}
					h_decoded_word[(time_count-ITERATE_TIME*4+1)%(ITERATE_TIME*4)][c] = decision;
				}
			}
		}
	}
};
#endif

#ifdef CODE2
extern "C"
void fun_matrix()
{
	// checkMatrix
	int row_deg[12] = { 3, 3, 3, 3, 10, 10, 10, 10, 10, 10, 10, 10 }; //一共12个行块

	check_node c_node[120] =
	{
		{ 9, 0, 0 }, { 17, 0, 0 }, { 20, 1, 160 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 10, 0, 0 }, { 17, 1, 0 }, { 18, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 11, 0, 0 }, { 18, 2, 0 }, { 19, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 5, 0, 0 }, { 12, 2, 0 }, { 13, 2, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 6, 2, 0 }, { 10, 2, 214 }, { 11, 3, 185 }, { 12, 3, 251 }, { 14, 0, 248 }, { 15, 0, 12 }, { 16, 0, 111 }, { 18, 0, 0 }, { 22, 0, 0 }, { 26, 0, 0 },
		{ 7, 1, 0 }, { 11, 4, 182 }, { 12, 4, 249 }, { 13, 3, 65 }, { 15, 1, 55 }, { 16, 1, 12 }, { 17, 0, 227 }, { 19, 0, 0 }, { 23, 0, 0 }, { 27, 0, 0 },
		{ 1, 2, 0 }, { 3, 3, 214 }, { 5, 5, 35 }, { 6, 4, 167 }, { 7, 1, 23 }, { 9, 2, 147 }, { 10, 1, 54 }, { 13, 0, 0 }, { 17, 0, 0 }, { 21, 0, 0 },
		{ 2, 2, 0 }, { 3, 4, 7 }, { 4, 5, 31 }, { 6, 5, 162 }, { 7, 2, 99 }, { 8, 2, 105 }, { 10, 2, 133 }, { 14, 0, 0 }, { 18, 0, 0 }, { 22, 0, 0 },
		{ 0, 2, 184 }, { 3, 5, 0 }, { 7, 3, 0 }, { 11, 1, 66 }, { 13, 1, 173 }, { 14, 1, 42 }, { 15, 1, 0 }, { 21, 1, 209 }, { 22, 1, 103 }, { 27, 0, 90 },
		{ 1, 3, 0 }, { 4, 2, 243 }, { 5, 1, 42 }, { 7, 2, 52 }, { 9, 1, 0 }, { 12, 1, 141 }, { 15, 2, 70 }, { 21, 0, 237 }, { 22, 0, 77 }, { 25, 0, 0 },
		{ 2, 3, 0 }, { 4, 3, 20 }, { 5, 2, 197 }, { 6, 2, 93 }, { 10, 1, 0 }, { 12, 2, 84 }, { 13, 1, 206 }, { 22, 1, 122 }, { 23, 0, 67 }, { 26, 0, 0 },
		{ 3, 3, 0 }, { 5, 3, 97 }, { 6, 3, 91 }, { 7, 3, 17 }, { 11, 1, 0 }, { 13, 2, 164 }, { 14, 2, 11 }, { 20, 1, 125 }, { 23, 1, 237 }, { 27, 0, 0 }
	};  //不够最大行重进行填充


	for (int i = 0; i<12; i++)
	{
		h_matrix_node_c[i].row_deg = row_deg[i];
		for (int j = 0; j<10; j++)
		{
			h_matrix_node_c[i].order[j] = c_node[i * 10 + j];
		}
	}

	//variableMatrix
	int col_deg[28] = { 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 3, 3, 3, 3, 6, 6, 6, 6 }; //一共20个列块

	variable_node v_node[168] =
	{
		{ 4, 4, 248 }, { 6, 4, 23 }, { 7, 4, 99 }, { 8, 2, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 4, 5, 12 }, { 5, 4, 55 }, { 7, 5, 105 }, { 9, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 4, 6, 111 }, { 5, 5, 12 }, { 6, 5, 147 }, { 10, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 5, 6, 227 }, { 6, 6, 54 }, { 7, 6, 133 }, { 11, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 4, 7, 0 }, { 8, 3, 66 }, { 9, 1, 243 }, { 10, 1, 20 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 5, 7, 0 }, { 9, 2, 42 }, { 10, 2, 197 }, { 11, 1, 97 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 6, 7, 0 }, { 8, 4, 173 }, { 10, 3, 93 }, { 11, 2, 91 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 4, 7, 0 }, { 5, 5, 42 }, { 6, 3, 52 }, { 8, 3, 17 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 1, 8, 0 }, { 5, 6, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 2, 8, 0 }, { 6, 4, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 3, 8, 0 }, { 7, 4, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 4, 8, 0 }, { 8, 4, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 1, 9, 0 }, { 6, 5, 141 }, { 7, 5, 84 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 2, 9, 0 }, { 7, 6, 206 }, { 8, 5, 164 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 0, 9, 0 }, { 2, 7, 209 }, { 5, 6, 11 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 1, 9, 0 }, { 2, 8, 103 }, { 3, 6, 70 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 6, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 7, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 8, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 9, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 2, 9, 90 }, { 5, 7, 125 }, { 10, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 0, 7, 237 }, { 8, 0, 0 }, { 11, 0, 184 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 0, 8, 77 }, { 1, 7, 122 }, { 9, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 1, 8, 67 }, { 2, 8, 237 }, { 10, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 },
		{ 3, 1, 0 }, { 4, 1, 0 }, { 7, 1, 241 }, { 9, 1, 214 }, { 10, 1, 7 }, { 11, 1, 0 },
		{ 0, 9, 0 }, { 4, 2, 0 }, { 5, 1, 0 }, { 7, 2, 185 }, { 8, 1, 182 }, { 10, 2, 31 },
		{ 1, 9, 0 }, { 5, 2, 0 }, { 6, 1, 0 }, { 7, 3, 251 }, { 8, 2, 249 }, { 9, 2, 35 },
		{ 2, 9, 0 }, { 3, 2, 160 }, { 6, 2, 0 }, { 8, 3, 65 }, { 9, 3, 167 }, { 10, 3, 162 }
	}; //不够最大列重进行填充


	for (int i = 0; i<28; i++)
	{
		h_matrix_node_v[i].col_deg = col_deg[i];
		for (int j = 0; j<6; j++)
		{
			h_matrix_node_v[i].order[j] = v_node[i * 6 + j];
		}
	}
}

extern "C"
void update(int time_count, float *h_channel_info, int **h_decoded_word, INFO_COL *h_info_col_2_row, INFO_ROW *h_info_row_2_col)
{
	int row_offset = time_count * 768;
	int col_offset = (time_count - 3) * 1792;
	int row_id = 0;
	int col_id = 0;
	node_position node[10];
	int relative_id = 0;
	int block_id = 0;
	int deg;
	int i = 0;
	int number = 0;
	int offset = 0;

	float info[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int info_symbol[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	float info_sum = 0;
	int info_sum_symbol = 0;
	float info_channel = 0;

	//Check Update

	for (int I = 0; I<ITERATE_TIME; I++)
	{
		for (int r = 0; r<CHECK_SIZE; r++)
		{
			row_id = row_offset - I * 4 * 768 + r;
			i = 0;
			info_sum = 0;
			info_sum_symbol = 0;

			if (row_id >= 0)
			{
				relative_id = row_id % 256;
				block_id = (row_id % 3072) / 256;
				deg = h_matrix_node_c[block_id].row_deg;
				offset = ((row_id / 768) - ROW_OFFSET) * 1792;
				while (i<deg)
				{
					node[i].pos = h_matrix_node_c[block_id].order[i].col;
					node[i].cpos = h_matrix_node_c[block_id].order[i].col_2_row;
					number = h_matrix_node_c[block_id].order[i].number;
					node[i].pos = offset + node[i].pos * 256 + (number + relative_id) % 256;
					node[i].pos = node[i].pos%COL_LENGTH;


					if (node[i].pos >= 0)
					{
						info[i] = h_info_col_2_row->info[node[i].cpos][node[i].pos];
					}
					else
					{
						info[i] = 0;
					}

					i++;
				}

				i = 0;
				while (i<deg)
				{
					unsigned int* pa = (unsigned int *)(&info[i]);
					info_symbol[i] = *pa >> 31;
					if (info[i] != 0)
					{
						info[i] = dev_ltanh(info[i]);
						info_sum += info[i];
						info_sum_symbol ^= info_symbol[i];
					}
					else
					{
						info[i] = 0;
					}
					i++;
				}

				i = 0;
				while (i<deg)
				{
					if (info[i] != 0)
					{
						info[i] = info_sum - info[i];
						info[i] = dev_ltanh(info[i]);
						unsigned int* pa = (unsigned int *)(&info[i]);
						*pa ^= ((info_sum_symbol^info_symbol[i]) << 31);

						h_info_row_2_col->info[i][row_id%ROW_LENGTH] = info[i];
					}
					i++;
				}

			}
		}
	}




	//Variable Update

	for (int I = 0; I<ITERATE_TIME; I++)
	{
		for (int c = 0; c<BLOCK_SIZE; c++)
		{
			col_id = col_offset - I * 4 * 1792 + c;
			i = 0;
			info_sum = 0;
			info_sum_symbol = 0;
			info_channel = 0;

			if (col_id >= 0)
			{
				relative_id = col_id % 256;
				block_id = (col_id % 7168) / 256;
				deg = h_matrix_node_v[block_id].col_deg;
				number = 0;
				offset = (col_id / 1792) * 768;

				info_channel = h_channel_info[col_id%COL_LENGTH];

				while (i<deg)
				{
					node[i].pos = h_matrix_node_v[block_id].order[i].row;
					node[i].cpos = h_matrix_node_v[block_id].order[i].row_2_col;
					number = h_matrix_node_v[block_id].order[i].number;
					node[i].pos = offset + node[i].pos * 256 + (relative_id - number + 256) % 256;
					node[i].pos = node[i].pos%ROW_LENGTH;

					info[i] = h_info_row_2_col->info[node[i].cpos][node[i].pos];

					i++;
				}

				i = 0;
				while (i<deg)
				{
					info_sum += info[i];
					i++;
				}

				i = 0;
				while (i<deg)
				{
					float test_info = info[i];
					float test_info_sum = info_sum;
					info[i] = info_sum - info[i] + info_channel;

					h_info_col_2_row->info[i][col_id%COL_LENGTH] = info[i];
					i++;
				}


				//Make decision

				if (I == (ITERATE_TIME - 1))
				{
					int decision = 0;

					info_sum += info_channel;
					if (info_sum<0)
					{
						decision = 1;
					}
					else
					{
						decision = 0;
					}
					h_decoded_word[(time_count - ITERATE_TIME * 4 + 1) % (ITERATE_TIME * 4)][c] = decision;
				}
			}
		}
	}
};
#endif