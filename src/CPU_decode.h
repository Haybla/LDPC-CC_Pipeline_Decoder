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

extern "C"
{
	void fun_matrix();
	int countError(int *word);
	float dev_ltanh(float info);
	void update(int time_count, float *h_channel_info, int **h_decoded_word, INFO_COL *h_info_col_2_row, INFO_ROW *h_info_row_2_col);
}