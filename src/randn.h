#define PI 3.141592653

double Normal(double x,double miu,double sigma) //概率密度函数
{
	return 1.0/sqrt(2*PI*sigma)*exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
}
double AverageRandom(double min,double max)
{
	double temp;
	temp=rand();
	temp = temp/(double)RAND_MAX;//0~1
	temp=temp*(max-min)+min;//min~max
	//temp=2*max*temp;//-max~+max
	return temp;
}
double NormalRandom(double miu,double sigma,double min,double max)//产生正态分布随机数
{
	double x;
	double dScope;
	double y;
	double P;

	do
	{
		x = AverageRandom(min,max);//产生min和max之间的随机数
		y = Normal(x, miu, sigma);//x点处的概率大小y
		P=Normal(miu,miu,sigma);
		dScope = AverageRandom(0,P);//x=0点的概率大小dScope
	}while( dScope > y);
	return x;
}
float randn(int size,float *rand_N)
{
	int i;
	//char s[10];
	//string s;
	//rand_N = (float*)malloc(sizeof(float)*size);

	for(i=0;i<size;i++)
	{
		rand_N[i]=(float)NormalRandom(0,1,-6,+6);
	}
	//把数据输出到文件中

	return 0;
}
double gaussrand()
{   
	double n=0;   
	for(int i=0;i<12;i++)   
	{   
		n+=(double)rand()/RAND_MAX;   
	}   
	n=(n-6); //标准化   
	return n;
} ;
float NoCal(float dB)
{
	float SNRpbit,No_uncoded,R,No;
	//SNRpbit= exp ((dB/10)*log(10.));
	SNRpbit = pow((float)10.0,(float)(dB/10.0));
	No_uncoded=(float)1.0/SNRpbit;
	R=(float)((float)MSG_SIZE/(float)BLOCK_SIZE);
	//R = 0.4;
	No=No_uncoded/R;
	return No;
};
void V_rand(float No,float *fV)
{
	double sigma;
	float temp_fv[BLOCK_SIZE];// 
	//temp_fv = (float*)malloc(sizeof(float)*BLOCK_SIZE);
	randn(BLOCK_SIZE,temp_fv);
	sigma=sqrt(No/2.0);
	for (int i=0;i<BLOCK_SIZE;i++)
	{
		//fV[i] =sigma*temp_fv[i];
		fV[i]= float(-1.0) + sigma*temp_fv[i];  //map 0->-1,1->+1
		fV[i] *= (((float)4.0/No)*(-1));
	}
}