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

#define PI 3.141592653

double Normal(double x,double miu,double sigma)
{
	return 1.0/sqrt(2*PI*sigma)*exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
};

double AverageRandom(double min,double max)
{
	double temp;
	temp=rand();
	temp = temp/(double)RAND_MAX;
	temp=temp*(max-min)+min;
	return temp;
};

double NormalRandom(double miu,double sigma,double min,double max)
{
	double x;
	double dScope;
	double y;
	double P;

	do
	{
		x = AverageRandom(min,max);
		y = Normal(x, miu, sigma);
		P=Normal(miu,miu,sigma);
		dScope = AverageRandom(0,P);
	}while( dScope > y);
	return x;
};

float randn(int size,float *rand_N)
{
	int i;

	for(i=0;i<size;i++)
	{
		rand_N[i]=(float)NormalRandom(0,1,-6,+6);
	}

	return 0;
};

double gaussrand()
{   
	double n=0;   
	for(int i=0;i<12;i++)   
	{   
		n+=(double)rand()/RAND_MAX;   
	}   
	n=(n-6);
	return n;
};

float NoCal(float dB)
{
	float SNRpbit,No_uncoded,R,No;
	SNRpbit = pow((float)10.0,(float)(dB/10.0));
	No_uncoded=(float)1.0/SNRpbit;
	R=(float)((float)MSG_SIZE/(float)BLOCK_SIZE);
	No=No_uncoded/R;
	return No;
};

void V_rand(float No,float *fV)
{
	double sigma;
	float temp_fv[BLOCK_SIZE];// 
	randn(BLOCK_SIZE,temp_fv);
	sigma=sqrt(No/2.0);
	for (int i=0;i<BLOCK_SIZE;i++)
	{
		fV[i]= float(1.0) + sigma*temp_fv[i];  //map 0->1,1->-1
		fV[i] *= (float)4.0/No;
	}
};