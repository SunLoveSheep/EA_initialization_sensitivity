#include "stdafx.h"
//#include "cec14_eotp.h"
//#include <iostream>
//#include <math.h>
#//include <malloc.h>
//#include <memory.h>

using namespace std;

cec14_eotp::cec14_eotp()
{
	
}

cec14_eotp::~cec14_eotp()
{
	
}

double *OShift_cec14eotp,*M_cec14eotp,*y_cec14eotp,*z_cec14eotp;
int ini_flag_cec14eotp=0,n_flag_cec14eotp,func_flag_cec14eotp;

//#define e  2.7182818284590452353602874713526625
//#define pi 3.1415926535897932384626433832795029

void cec14_eotp::FileIO(int dim, int func_num)
{
	int cf_num=10,i;
	if (ini_flag_cec14eotp==1)
	{
		if ((n_flag_cec14eotp!=dim)||(func_flag_cec14eotp!=func_num))
		{
			ini_flag_cec14eotp=0;
		}
	}

	if (ini_flag_cec14eotp==0)
	{
		FILE *fpt;
		char FileName[30];
		
		y_cec14eotp=(double *)malloc(sizeof(double)  *  dim);//dim: problem dimension
		z_cec14eotp=(double *)malloc(sizeof(double)  *  dim);
		
		int func_num_read = (int) ceil(func_num / 3.0);

		//sprintf(FileName, "input_data/M_%d_D%d.txt",func_num_read,dim);
		sprintf_s(FileName, "input_data_cec14_eotp/M_%d_D%d.txt",func_num_read,dim);
		//fpt = fopen(FileName,"r");
		fopen_s(&fpt,FileName,"r");
		if (fpt==NULL)
		{
		    printf("\n Error: Cannot open input file for reading \n");
		}

		M_cec14eotp=(double*)malloc(cf_num*dim*dim*sizeof(double));
		if (M_cec14eotp==NULL)
			printf("\nError: there is insufficient memory available!\n");
		for (i=0; i<cf_num*dim*dim; i++)
		{
			//fscanf(fpt,"%Lf",&M[i]);
			fscanf_s(fpt,"%Lf",&M_cec14eotp[i]);
			//printf("M[%d] = %LE,",i+1,M[i]);
		}
		fclose(fpt);
		
		
		//sprintf(FileName, "input_data/shift_data_%d.txt",func_num_read);
		sprintf_s(FileName, "input_data_cec14_eotp/shift_data_%d.txt",func_num_read);
		//fpt = fopen(FileName,"r");
		fopen_s(&fpt,FileName,"r");
		if (fpt==NULL)
		{
			printf("\n Error: Cannot open input file for reading \n");
		}
		OShift_cec14eotp=(double *)malloc(dim*cf_num*sizeof(double));
		if (OShift_cec14eotp==NULL)
			printf("\nError: there is insufficient memory available!\n");
		for(i=0;i<cf_num*dim;i++)
		{
			//fscanf(fpt,"%Lf",&OShift[i]);
			fscanf_s(fpt,"%Lf",&OShift_cec14eotp[i]);
		}
		fclose(fpt);

		n_flag_cec14eotp=dim;
		func_flag_cec14eotp=func_num;
		ini_flag_cec14eotp=1;
		printf("Function has been initialized!\n");
	}
}

void cec14_eotp::Release()
{
	free(M_cec14eotp);//Rotation matrix, size n*n
	free(OShift_cec14eotp);//shifting vector, size n
	free(y_cec14eotp);//shifting resultant vector
	//free(z);//rotation (if has) resultant
}

/* 1 Sphere function*/
double cec14_eotp::sphere(const double *x, int dim) {
	int i;
	double ret;

	ret = 0.0;
	for (i = 0; i < dim; i++) {
		ret += (x[i] * x[i]);
	}
	//cout<<ret<<endl;
	//system("pause");
	return ret;
}

/* 2 Ellipsoid function*/
double cec14_eotp::ellipsoid(const double *x, int dim) {
	int i;
	double ret;

	ret = 0.0L;
	for (i = 0; i < dim; i++) {
		ret += (x[i] * x[i] * (i + 1.0L));
	}

	return ret;
}

/* 3 Rotated Ellipsoid function*/
double cec14_eotp::rotated_ellipsoid_revised(const double *x, int dim) {
	double * y;
	int i, j;
	double ret;

	y = (double*) malloc(sizeof(double) * dim);
	memset(y, 0, dim * sizeof(double));
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			y[i] += M_cec14eotp[i * dim + j] * x[j];
		}
	}

	ret = ellipsoid(y, dim);

	free(y);
	/*if (dim != 2 && dim != 10 && dim != 20 && dim != 30)
		free(M);*/
	return ret;
}

/* 4 Step function*/
double cec14_eotp::step(const double *x, int dim) {
	int i;
	double tmp, ret;

	ret = 0.0L;
	for (i = 0; i < dim; i++) {
		tmp = floor(x[i] + 0.5L);
		ret += (tmp * tmp);
	}

	return ret;
}

/* 5 Ackley function*/
double cec14_eotp::ackley(const double *x, int dim) {
	int i;
	double sum1, sum2;
	double ret;

	sum1 = 0.0L;
	sum2 = 0.0L;
	for (i = 0; i < dim; i++) {
		sum1 += (x[i] * x[i]);
		sum2 += cos(2.0L * pi * x[i]);
	}

	ret = -20.0L * exp(-0.2L * sqrt(sum1 / dim)) - exp(sum2 / dim) + 20.0L + e;

	return ret;
}

/* 6 Griewank function*/
double cec14_eotp::griewank(const double *x, int dim) {
	double sum, prod;
	double ret;
	int i;

	sum = 0.0L;
	prod = 1.0L;

	for (i = 0; i < dim; i++) {
		sum += x[i] * x[i];
		prod *= cos(x[i] / sqrt(i + 1.0L));
	}

	ret = 1.0L + sum / 4000.0L - prod;

	return ret;
}

/* 7 Rosenbrock function*/
double cec14_eotp::rosenbrock(const double *x, int dim) {
	double tmp1, tmp2;
	double ret;
	int i;

	ret = 0.0L;
	for (i = 0; i < dim - 1; i++) {
		tmp1 = x[i] * x[i] - x[i + 1];
		tmp2 = x[i] - 1.0L;
		ret += (100.0L * tmp1 * tmp1 + tmp2 * tmp2);
	}

	return ret;
}

/* 8 Rastrigin function*/
double cec14_eotp::rastrigin(const double *x, int dim) {
	int i;
	double ret;

	ret = 0.0L;
	for (i = 0; i < dim; i++) {
		ret += (x[i] * x[i] - 10.0L * cos(2.0L * pi * x[i]) + 10.0L);
	}
	return ret;
}

void cec14_eotp::scale(const double *x, double s, double * y, int n) {
	for (int i = 0; i < n; i++)
		y[i] = x[i] * s;
}

/*inline void shift(const double* x,const double * shift, double * y, int n) {
	for (int i = 0; i < n; i++)
		y[i] = shift[i] + x[i];
}*/

void cec14_eotp::shiftfunc (const double *x, double *xshift, int dim,const double *Os)
{
	int i;
    for (i=0; i<dim; i++)
    {
        xshift[i]=x[i]-Os[i];
    }
}

void cec14_eotp::shift_s(const double* x, double shift, double* y, int n) {
	for (int i = 0; i < n; i++)
		y[i] = shift + x[i];
}

void cec14_eotp::m_mul(const double * M, const double * x, double *y, int n) {
	double *z;
	z = (double*) malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) {
		z[i] = 0;
		for (int j = 0; j < n; j++) {
			z[i] += M[i * n + j] * x[j];
		}
	}

	for (int i = 0; i < n; i++)
		y[i] = z[i];

	free(z);
}

double cec14_eotp::cec14_eotp_problems(const double * x, int dim, int no) 
{
	int func_no;
	int dim_supposed;
	double * y;
	//const double * o = new double[dim];
	//const double * M = new double[dim*dim];
	double ret=0;
	//bool isGood = true;

	if (no < 1 || no > 24) {
		//isGood = false;
		cout<<"Wrong Function Number! Should be 1-24";
		return 0.0;
	}

	switch (no % 3) {
	case 1:
		dim_supposed = 10;
		break;
	case 2:
		dim_supposed = 30;
		break;
	case 0:
		dim_supposed = 100;
		break;
	}

	if (dim != dim_supposed) {
		//*isGood = false;
		cout<<"Wrong Dimension! Should be 10/30/100";
		getchar();
		return 0.0;
	}

	func_no = (int) ceil(no / 3.0);
	//M = getMatrix(func_no, dim);
	//o = getO(func_no);
	
	y = (double*) malloc(dim * sizeof(double));

	switch (func_no) {
	case 1:
		scale(OShift_cec14eotp, 1.0 / 10.0, y, dim);
		//shift(y, x, y, dim);
		shiftfunc(x, y, dim, y);
		ret = sphere(y, dim);
		break;
	case 2:
		scale(OShift_cec14eotp, 1.0 / 10.0, y, dim);
		shiftfunc(x, y, dim, y);

		ret = ellipsoid(y, dim);
		break;

	case 3:
		scale(OShift_cec14eotp, 1.0 / 10.0, y, dim);
		shiftfunc(x, y, dim, y);
		//m_mul(M, y, y, dim);
		//ret = rotated_ellipsoid(y, dim);

		ret = rotated_ellipsoid_revised(y, dim);
		break;

	case 4:
		scale(OShift_cec14eotp, 1.0 / 10.0, y, dim);
		shiftfunc(x, y, dim, y);

		ret = step(y, dim);
		break;

	case 5:
		scale(OShift_cec14eotp, 1.0 / 10.0, y, dim);
		shiftfunc(x, y, dim, y);

		ret = ackley(y, dim);
		break;
	case 6:
		scale(OShift_cec14eotp, 1.0 / 10.0, y, dim);
		shiftfunc(x, y, dim, y);

		ret = griewank(y, dim);
		break;

	case 7:
		scale(OShift_cec14eotp, 1.0 / 10.0, y, dim);
		shiftfunc(x, y, dim, y);
		m_mul(M_cec14eotp, y, y, dim);
		scale(y, 2.08 / 20, y, dim);
		shift_s(y, 1, y, dim);

		ret = rosenbrock(y, dim);
		break;

	case 8:
		scale(OShift_cec14eotp, 1.0 / 10.0, y, dim);
		shiftfunc(x, y, dim, y);
		m_mul(M_cec14eotp, y, y, dim);
		scale(y, 5.12 / 20, y, dim);
		ret = rastrigin(y, dim);
		break;
	}

	extern Solution solution;
	for (int i=0;i<dim;i++)
	{
		solution.GlobalOptima[i] = OShift_cec14eotp[i];
	}

	free(y);//shifting resultant vector
	return ret;
}