#pragma once

#include "totalDefine.h"

extern "C"
{
	void fun_matrix();
	int countError(int *word);
	float dev_ltanh(float info);
	void update(int time_count, float *h_channel_info, int **h_decoded_word, INFO_COL *h_info_col_2_row, INFO_ROW *h_info_row_2_col);
}