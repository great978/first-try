// commonFFT.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。

#include "pch.h"
#define _CRT_SECURE_NO_WARNINGS
#include "stdlib.h"
#include "math.h"
#include "stdio.h"
#include "time.h"
#include "Windows.h"

#define Pi  3.1416   //定义圆周率Pi
#define LEN sizeof(struct Compx)  //定义复数结构体大小

int readdata();
void FFT(struct Compx *Xin, struct Compx *Xout, int p[], int l, int HH[], int KK[]);
void DD(struct Compx *Xin, struct Compx *Xout, int p[], int HH[], int KK[], int l, struct Compx ** Val, int m);
void IFFT(struct Compx *Xin, struct Compx *Xout, int p[], int l, int HH[], int KK[]);
void IDD(struct Compx *Xin, struct Compx *Xout, int p[], int HH[], int KK[], int l, struct Compx ** Val, int m);

int data[1000];
int N = 0;
int m, count = 0;

//-----定义复数结构体-----------------------//
struct Compx
{
	float real;
	float imag;
} Compx;
//-----复数结构体清零函数-------------------//
void ZeroCompx(struct Compx X[], unsigned int KK)
{
	unsigned int i;
	for (i = 0; i < KK; i++)
	{
		X[i].imag = 0;
		X[i].real = 0;

	}
	return;
}
//-----复数乘法运算函数---------------------//
struct Compx EE(struct Compx b1, struct Compx b2)
{
	struct Compx b3;
	b3.real = b1.real*b2.real - b1.imag*b2.imag;
	b3.imag = b1.real*b2.imag + b1.imag*b2.real;
	return(b3);
}
struct Compx EE1(struct Compx b1, float b2)
{
	struct Compx b3;
	b3.real = b1.real*b2;
	b3.imag = b1.imag*b2;
	return(b3);
}
//-----复数加法运算函数---------------------//
struct Compx add(struct Compx a, struct Compx b)
{
	struct Compx c;
	c.real = a.real + b.real;
	c.imag = a.imag + b.imag;
	return(c);
}

int main()
{	/*-----------------------开启计时器--------------------------------*/

	/*-----------------------读入数据--------------------------------*/
	int i, j;
	N = readdata();
	/*-----------------------素数分解--------------------------------*/
	printf("长度为%d，素数分解：", N);
	m = 0;
	int n = N;
	int a[1000];                               //分解后的数组
	printf("%d=", n);
	for (i = 2; n >= i; i++)
	{
		while (n%i == 0)
		{
			printf("%d", i);
			a[m] = i;
			n /= i;
			if (n >= 1)
			{
				printf("*");
				m++;
			}
		}
	}
	if (n >= 1) {
		printf("%d", n);
		a[m] = n;
	}
	printf("\n");
	int * p = (int *)malloc(m * sizeof(int));
	m += 1;
	/*-----------------------排序--------------------------------*/
	for (i = N - 1; i > 0; --i)          //从n-1循环的到0，也是n次
		for (j = 0; j < i; j++)
			if (a[j] < a[j + 1])
			{
				int temp = a[j];
				a[j] = a[j + 1];
				a[j + 1] = temp;
			}
	for (i = 0; i < m; i++) {
		p[i] = a[i];
	}
	for (i = 0; i < m; i++) {
		printf("%d ", p[i]);
	}
	printf("\n");
	/*------------------------FFT变换-------------------------------*/
	/*
	struct   Compx Source[N];  //定义FFT的采样点存放数组
	struct   Compx Result[N];  //定义FFT的运算结果存放数组
	struct   Compx IResult[N];
	*/
	struct   Compx * Source = (struct Compx *)malloc(N*LEN);			//为结构体分配存储空间
	struct   Compx * Result = (struct Compx *)malloc(N*LEN);
	struct   Compx * IResult = (struct Compx *)malloc(N*LEN);
	/*int HH[m];
	int KK[m];
	int MM[m];
	int YY[m];*/
	int * HH = (int *)malloc(m * sizeof(int));    //为数组分配内存空间
	int * KK = (int *)malloc(m * sizeof(int));
	int * MM = (int *)malloc(m * sizeof(int));
	int * YY = (int *)malloc(m * sizeof(int));

	printf("\nSource初始化：\n");
	for (i = 0; i < N; i++)
	{
		Source[i].real = data[i];
		Source[i].imag = 0;
		printf("%.4f ", Source[i].real);
		printf("+%.4fj  ", Source[i].imag);
		printf("\n");
	}
	printf("\Result初始化：\n");
	for (i = 0; i < N; i++)
	{
		Result[i].real = 0;
		Result[i].imag = 0;
		printf("%.4f ", Result[i].real);
		printf("+%.4fj  ", Result[i].imag);
		printf("\n");
	}
	printf("\IResult初始化：\n");
	for (i = 0; i < N; i++)
	{
		IResult[i].real = 0;
		IResult[i].imag = 0;
		printf("%.4f ", IResult[i].real);
		printf("+%.4fj  ", IResult[i].imag);
		printf("\n");
	}
	FFT(Source, Result, p, m, HH, KK);      //采用任意基数的FFT变换
	printf("\nFFT分析结果:\n");
	for (int i = 0; i < N; i++)
	{
		printf("%.4f", Result[i].real);
		printf("+%.4fj    ", Result[i].imag);
		printf("\n");
	}

	IFFT(Result, IResult, p, m, HH, KK);

	printf("\n");
	printf("\nIFFT分析结果:\n");
	for (int i = 0; i < N; i++)
	{
		printf("%.4f", IResult[i].real);
		printf("+%.4fj    ", IResult[i].imag);
		printf("\n");
	}

	printf("\n运行结束n");
	system("pause");
	system("pause");
	return 0;
}

void IFFT(struct Compx *Xin, struct Compx *Xout, int p[], int l, int HH[], int KK[])
{
	if (l > 1)
	{
		int *j = (int *)malloc((l - 1) * sizeof(int));
		for (j[l - 1] = 0; j[l - 1] < p[l - 1]; j[l - 1]++)
		{
			HH[l - 1] = j[l - 1];
			IFFT(Xin, Xout, p, l - 1, HH, KK);

		}
	}
	else
	{
		int *j = (int *)malloc(sizeof(int));
		for (j[0] = 0; j[0] < p[0]; j[0]++)
		{
			//struct Compx Val[m - 1][m];
			struct Compx ** Val;
			Val = (struct Compx**)malloc((m - 1)*LEN);
			for (int i = 0; i < (m - 1); i++) {

				Val[i] = (struct Compx*)malloc(m*LEN);

			}
			ZeroCompx(Val[0], m);
			HH[0] = j[0];
			DD(Xin, Xout, p, HH, KK, 0, Val, m);

			int jm = 0;
			for (int i = 1; i < m; i++)
			{
				int jjm = 1;
				for (int j = 0; j < i; j++)
				{
					jjm = jjm * p[j];
				}
				jm = jm + HH[i] * jjm;
			}
			jm = jm + HH[0];
			Xout[jm] = EE1(Xout[jm], 1.0 / N);
		}
	}
	if (l == m)
	{
		struct   Compx temp1 = Xout[0];
		for (int i = 0; i < N / 2; i++) {
			struct   Compx temp = Xout[i];
			Xout[i] = Xout[N - 1 - i];
			Xout[N - 1 - i] = temp;
		}
		for (int i = N - 1; i > 0; i--)
		{
			struct   Compx temp2 = Xout[i];
			Xout[i] = Xout[i - 1];
		}
		Xout[0] = temp1;

	}
}

void IDD(struct Compx *Xin, struct Compx *Xout, int p[], int HH[], int KK[], int l, struct Compx ** Val, int m)
{
	double ps, ps1;
	struct Compx  W, Temp;

	if (l < m - 1)
	{
		if (l > 0)
		{
			ZeroCompx(Val[l], m);
		}
		int * k = (int *)malloc(l * sizeof(int));
		for (k[l] = 0; k[l] < p[m - l - 1]; k[l]++)
		{
			KK[l] = k[l];
			l = l + 1;
			DD(Xin, Xout, p, HH, KK, l, Val, m);

			l = l - 1;
			int pad = 1;
			for (int i = m - l; i < m; i++)
			{
				pad = pad * p[i];
			}

			int jadd = 0;
			for (int i = 1; i <= m - 1 - l; i++)
			{
				int pnn = 1;
				for (int j = 0; j < i; j++)
				{
					pnn = pnn * p[j];
				}
				jadd = jadd + HH[i] * pnn;
			}
			jadd = jadd + HH[0];

			if (l == 0)
			{
				ps = jadd * k[l];
			}
			else
			{
				ps = jadd * k[l] * pad;
			}

			ps1 = -2 * Pi / N * ps;
			W.real = cos(ps1);
			W.imag = -sin(ps1);
			ps = 0;
			ps1 = 0;
			if (l >= 1)
			{
				Val[l][KK[l]] = EE(Val[l][KK[l]], W);

				Val[l - 1][KK[l - 1]] = add(Val[l][KK[l]], Val[l - 1][KK[l - 1]]);

			}
			else
			{
				Val[l][KK[l]] = EE(Val[l][KK[l]], W);
			}

			W.real = 0;
			W.imag = 0;
			int jm = 0;
			for (int i = 1; i < m; i++)
			{
				int jjm = 1;
				for (int j = 0; j < i; j++)
				{
					jjm = jjm * p[j];
				}
				jm = jm + HH[i] * jjm;
			}
			jm = jm + HH[0];

			if (l == 0)
			{
				Xout[jm] = add(Xout[jm], Val[l][KK[l]]);

			}
		}
	}
	else
	{
		int * k = (int *)malloc((m - 1) * sizeof(int));
		for (k[m - 1] = 0; k[m - 1] < p[0]; k[m - 1]++)
		{
			KK[m - 1] = k[m - 1];
			int padd = 1;
			for (int i = 1; i < m; i++)
			{
				padd = padd * p[i];
			}
			ps = padd * k[m - 1] * HH[0];
			ps1 = -2 * Pi / N * ps;
			W.real = cos(ps1);
			W.imag = -sin(ps1);
			int qm = 0;
			for (int i = 1; i < m; i++)
			{
				int pm = 1;
				for (int j = 0; j < i; j++)
				{
					pm = pm * p[m - j - 1];
				}
				qm = qm + KK[i] * pm;
			}
			qm = qm + KK[0];
			Temp = EE(Xin[qm], W);
			ps = 0;
			ps1 = 0;
			W.real = 0;
			W.imag = 0;
			Val[m - 2][KK[m - 2]] = add(Val[m - 2][KK[m - 2]], Temp);

		}
	}
}

void FFT(struct Compx *Xin, struct Compx *Xout, int p[], int l, int HH[], int KK[])
{
	if (l > 1)
	{
		int *j = (int *)malloc((l - 1) * sizeof(int));
		for (j[l - 1] = 0; j[l - 1] < p[l - 1]; j[l - 1]++)
		{
			HH[l - 1] = j[l - 1];
			FFT(Xin, Xout, p, l - 1, HH, KK);
		}
	}
	else
	{
		int *j = (int *)malloc(sizeof(int));
		for (j[0] = 0; j[0] < p[0]; j[0]++)
		{
			//struct Compx Val[m - 1][m];				//分配Val[m - 1][m]内存
			struct Compx ** Val;
			Val = (struct Compx**)malloc((m - 1)*LEN);
			for (int i = 0; i < (m - 1); i++) {

				Val[i] = (struct Compx*)malloc(m*LEN);

			}
			ZeroCompx(Val[0], m);
			HH[0] = j[0];
			DD(Xin, Xout, p, HH, KK, 0, Val, m);
		}
	}
}
void DD(struct Compx *Xin, struct Compx *Xout, int p[], int HH[], int KK[], int l, struct Compx ** Val, int m)
{
	float ps, ps1;
	struct Compx  W, Temp;

	if (l < m - 1)
	{
		if (l > 0)
		{
			ZeroCompx(Val[l], m);
		}
		int * k = (int *)malloc(l * sizeof(int));
		for (k[l] = 0; k[l] < p[m - l - 1]; k[l]++)
		{
			KK[l] = k[l];
			l = l + 1;
			DD(Xin, Xout, p, HH, KK, l, Val, m);

			l = l - 1;
			int pad = 1;
			for (int i = m - l; i < m; i++)
			{
				pad = pad * p[i];
			}

			int jadd = 0;
			for (int i = 1; i <= m - 1 - l; i++)
			{
				int pnn = 1;
				for (int j = 0; j < i; j++)
				{
					pnn = pnn * p[j];
				}
				jadd = jadd + HH[i] * pnn;
			}
			jadd = jadd + HH[0];

			if (l == 0)
			{
				ps = jadd * k[l];
			}
			else
			{
				ps = jadd * k[l] * pad;
			}

			ps1 = 2 * Pi / N * ps;
			W.real = cos(ps1);
			W.imag = -sin(ps1);
			ps = 0;
			ps1 = 0;
			if (l >= 1)
			{
				Val[l][KK[l]] = EE(Val[l][KK[l]], W);

				Val[l - 1][KK[l - 1]] = add(Val[l][KK[l]], Val[l - 1][KK[l - 1]]);

			}
			else
			{
				Val[l][KK[l]] = EE(Val[l][KK[l]], W);
			}

			W.real = 0;
			W.imag = 0;
			int jm = 0;
			for (int i = 1; i < m; i++)
			{
				int jjm = 1;
				for (int j = 0; j < i; j++)
				{
					jjm = jjm * p[j];
				}
				jm = jm + HH[i] * jjm;
			}
			jm = jm + HH[0];

			if (l == 0)
			{
				Xout[jm] = add(Xout[jm], Val[l][KK[l]]);

			}
		}
	}
	else
	{
		int * k = (int *)malloc((m - 1) * sizeof(int));
		for (k[m - 1] = 0; k[m - 1] < p[0]; k[m - 1]++)
		{
			KK[m - 1] = k[m - 1];
			int padd = 1;
			for (int i = 1; i < m; i++)
			{
				padd = padd * p[i];
			}
			ps = padd * k[m - 1] * HH[0];
			ps1 = 2 * Pi / N * ps;
			W.real = cos(ps1);
			W.imag = -sin(ps1);
			int qm = 0;
			for (int i = 1; i < m; i++)
			{
				int pm = 1;
				for (int j = 0; j < i; j++)
				{
					pm = pm * p[m - j - 1];
				}
				qm = qm + KK[i] * pm;
			}
			qm = qm + KK[0];
			Temp = EE(Xin[qm], W);
			ps = 0;
			ps1 = 0;
			W.real = 0;
			W.imag = 0;
			Val[m - 2][KK[m - 2]] = add(Val[m - 2][KK[m - 2]], Temp);

		}
	}
}

int readdata() {
	char c;
	int i = 0;
	int n = 0;
	while ((c = getchar()) != '\n')
	{
		if (isdigit(c))
		{
			ungetc(c, stdin);       //将c送回输入流
			scanf("%d", &data[n++]);
		}
	}
	return n;			//返回数组长度n
}

