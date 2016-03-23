#include "stdafx.h"
//#include "cec13.h"
//#include <iostream>
//#include <math.h>
//#include <malloc.h>
//#include <memory.h>

using namespace std;

cec13::cec13()
{
	
}

cec13::~cec13()
{
	
}

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag;

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

void cec13::FileIO(int dim, int func_num)
{
	int cf_num=10,i;
	if (ini_flag==1)
	{
		if ((n_flag!=dim)||(func_flag!=func_num))
		{
			ini_flag=0;
		}
	}

	if (ini_flag==0)
	{
		FILE *fpt;
		char FileName[30];
		free(M);
		free(OShift);
		free(y);
		free(z);
		free(x_bound);
		y=(double *)malloc(sizeof(double)  * dim);
		z=(double *)malloc(sizeof(double)  * dim);
		x_bound=(double *)malloc(sizeof(double)  * dim);
		for (i=0; i<dim; i++)
			x_bound[i]=100.0;
		
		sprintf_s(FileName, "input_data_cec13/M_D%d.txt", dim);
		fopen_s(&fpt,FileName,"r");
		if (fpt==NULL)
		{
		    printf("\n Error: Cannot open input file for reading \n");
		}

		M=(double*)malloc(cf_num*dim*dim*sizeof(double));
		if (M==NULL)
			printf("\nError: there is insufficient memory available!\n");
		for (i=0; i<cf_num*dim*dim; i++)
		{
				fscanf_s(fpt,"%Lf",&M[i]);
				//printf("M[%d] = %LE,",i+1,M[i]);
		}
		fclose(fpt);
		

		fopen_s(&fpt,"input_data_cec13/shift_data.txt","r");
		if (fpt==NULL)
		{
			printf("\n Error: Cannot open input file for reading \n");
		}
		OShift=(double *)malloc(dim*cf_num*sizeof(double));
		if (OShift==NULL)
			printf("\nError: there is insufficient memory available!\n");
		for(i=0;i<cf_num*dim;i++)
		{
				fscanf_s(fpt,"%Lf",&OShift[i]);
		}
		fclose(fpt);

		n_flag=dim;
		func_flag=func_num;
		ini_flag=1;
		//printf("Function has been initialized!\n");
	}
}

void cec13::Release()
{
	free(M);//Rotation matrix, size n*n
	free(OShift);//shifting vector, size n
	free(y);//shifting resultant vector
	//free(z);//rotation (if has) resultant
}

void shiftfunc (const double *x, double *xshift, int nx,double *Os)
{
	int i;
    for (i=0; i<nx; i++)
    {
        xshift[i]=x[i]-Os[i];
    }
}

void rotatefunc (const double *x, double *xrot, int nx,double *Mr)
{
	int i,j;
    for (i=0; i<nx; i++)
    {
        xrot[i]=0;
			for (j=0; j<nx; j++)
			{
				xrot[i]=xrot[i]+x[j]*Mr[i*nx+j];
			}
    }
}

void asyfunc (const double *x, double *xasy, int nx, double beta)
{
	int i;
    for (i=0; i<nx; i++)
    {
		if (x[i]>0)
        xasy[i]=pow(x[i],1.0+beta*i/(nx-1)*pow(x[i],0.5));
    }
}

void oszfunc (const double *x, double *xosz, int nx)
{
	int i,sx;
	double c1,c2,xx;
    for (i=0; i<nx; i++)
    {
		if (i==0||i==nx-1)
        {
			if (x[i]!=0)
				xx=log(fabs(x[i]));
			if (x[i]>0)
			{	
				c1=10;
				c2=7.9;
			}
			else
			{
				c1=5.5;
				c2=3.1;
			}	
			if (x[i]>0)
				sx=1;
			else if (x[i]==0)
				sx=0;
			else
				sx=-1;
			xosz[i]=sx*exp(xx+0.049*(sin(c1*xx)+sin(c2*xx)));
		}
		else
			xosz[i]=x[i];
    }
}

double cf_cal(const double *x, int nx, double *Os,double * delta,double * bias,double * fit, int cf_num)
{
	int i,j;
	double *w;
	double w_max=0,w_sum=0;
	w=(double *)malloc(cf_num * sizeof(double));
	for (i=0; i<cf_num; i++)
	{
		fit[i]+=bias[i];
		w[i]=0;
		for (j=0; j<nx; j++)
		{
			w[i]+=pow(x[j]-Os[i*nx+j],2.0);
		}
		if (w[i]!=0)
			w[i]=pow(1.0/w[i],0.5)*exp(-w[i]/2.0/nx/pow(delta[i],2.0));
		else
			w[i]=INF;
		if (w[i]>w_max)
			w_max=w[i];
	}

	for (i=0; i<cf_num; i++)
	{
		w_sum=w_sum+w[i];
	}
	if(w_max==0)
	{
		for (i=0; i<cf_num; i++)
			w[i]=1;
		w_sum=cf_num;
	}
	double ret = 0.0;
    for (i=0; i<cf_num; i++)
    {
		ret=ret+w[i]/w_sum*fit[i];
    }
	free(w);

	return ret;
}

double sphere_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Sphere */
{
	double ret=0.0;
	int i;
	shiftfunc(x, y, nx, Os);
	/*for (int i=0;i<nx;i++)
	{
		cout<<Os[i]<<" ";
	}
	cout<<endl;
	getchar();*/
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];
    for (i=0; i<nx; i++)
    {
        ret += z[i]*z[i];
    }

	return ret;
}

double ellips_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Ellipsoidal */
{
    int i;
	shiftfunc(x, y, nx, Os);
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];
    oszfunc (z, y, nx);
	double ret = 0.0;
    for (i=0; i<nx; i++)
    {
        ret += pow(10.0,6.0*i/(nx-1))*y[i]*y[i];
    }

	return ret;
}

double bent_cigar_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Bent_Cigar */
{
    int i;
	double beta=0.5;
	shiftfunc(x, y, nx, Os);
	if (r_flag==1)
		rotatefunc(y, z, nx, Mr);
	else
		for (i=0; i<nx; i++)
			z[i]=y[i];

    asyfunc (z, y, nx,beta);
	if (r_flag==1)
		rotatefunc(y, z, nx, &Mr[nx*nx]);
	else
		for (i=0; i<nx; i++)
			z[i]=y[i];

	double ret = 0.0;
	ret = z[0]*z[0];
    for (i=1; i<nx; i++)
    {
        ret += pow(10.0,6.0)*z[i]*z[i];
    }

	return ret;
}

double discus_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Discus */
{
    int i;
	shiftfunc(x, y, nx, Os);
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];
    oszfunc (z, y, nx);

	double ret = pow(10.0,6.0)*y[0]*y[0];
    for (i=1; i<nx; i++)
    {
        ret += y[i]*y[i];
    }

	return ret;
}

double dif_powers_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Different Powers */
{
	int i;
	shiftfunc(x, y, nx, Os);
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];
	double ret = 0.0;
    for (i=0; i<nx; i++)
    {
        ret += pow(fabs(z[i]),2+4*i/(nx-1));
    }
	ret=pow(ret,0.5);

	return ret;
}

double rosenbrock_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Rosenbrock's */
{
    int i;
	double tmp1,tmp2;
	shiftfunc(x, y, nx, Os);//shift
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]=y[i]*2.048/100;
    }
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);//rotate
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];
	for (i=0; i<nx; i++)//shift to orgin
    {
        z[i]=z[i]+1;
    }

    double ret = 0.0;
    for (i=0; i<nx-1; i++)
    {
		tmp1=z[i]*z[i]-z[i+1];
		tmp2=z[i]-1.0;
        ret += 100.0*tmp1*tmp1 +tmp2*tmp2;
    }

	return ret;
}

double schaffer_F7_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Schwefel's 1.2  */
{
    int i;
	double tmp;
	shiftfunc(x, y, nx, Os);
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];
	asyfunc (z, y, nx, 0.5);
	for (i=0; i<nx; i++)
		z[i] = y[i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(z, y, nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[i]=z[i];

	for (i=0; i<nx-1; i++)
		z[i]=pow(y[i]*y[i]+y[i+1]*y[i+1],0.5);
    double ret = 0.0;
    for (i=0; i<nx-1; i++)
    {
	  tmp=sin(50.0*pow(z[i],0.2));
      ret += pow(z[i],0.5)+pow(z[i],0.5)*tmp*tmp ;
    }
	ret = ret*ret/(nx-1)/(nx-1);
	
	return ret;
}

double ackley_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Ackley's  */
{
    int i;
    double sum1, sum2;

	shiftfunc(x, y, nx, Os);
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

	asyfunc (z, y, nx, 0.5);
	for (i=0; i<nx; i++)
		z[i] = y[i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(z, y, nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[i]=z[i];

    sum1 = 0.0;
    sum2 = 0.0;
    for (i=0; i<nx; i++)
    {
        sum1 += y[i]*y[i];
        sum2 += cos(2.0*PI*y[i]);
    }
    sum1 = -0.2*sqrt(sum1/nx);
    sum2 /= nx;
    double ret =  E - 20.0*exp(sum1) - exp(sum2) +20.0;

	return ret;
}

double weierstrass_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Weierstrass's  */
{
    int i,j,k_max;
    double sum,sum2, a, b;

	shiftfunc(x, y, nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]=y[i]*0.5/100;
    }
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

	asyfunc (z, y, nx, 0.5);
	for (i=0; i<nx; i++)
		z[i] = y[i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(z, y, nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[i]=z[i];

    a = 0.5;
    b = 3.0;
    k_max = 20;
    double ret = 0.0;
    for (i=0; i<nx; i++)
    {
        sum = 0.0;
		sum2 = 0.0;
        for (j=0; j<=k_max; j++)
        {
            sum += pow(a,j)*cos(2.0*PI*pow(b,j)*(y[i]+0.5));
			sum2 += pow(a,j)*cos(2.0*PI*pow(b,j)*0.5);
        }
        ret += sum;
    }
	ret -= nx*sum2;

	return ret;
}

double griewank_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Griewank's  */
{
    int i;
    double s, p;

	shiftfunc(x, y, nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]=y[i]*600.0/100.0;
    }
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

	for (i=0; i<nx; i++)
		z[i] = z[i]*pow(100.0,1.0*i/(nx-1)/2.0);


    s = 0.0;
    p = 1.0;
    for (i=0; i<nx; i++)
    {
        s += z[i]*z[i];
        p *= cos(z[i]/sqrt(1.0+i));
    }
    double ret = 1.0 + s/4000.0 - p;

	return ret;
}

double rastrigin_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Rastrigin's  */
{
    int i;
	double alpha=10.0,beta=0.2;
	shiftfunc(x, y, nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]=y[i]*5.12/100;
    }

	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

    oszfunc (z, y, nx);
    asyfunc (y, z, nx, beta);

	if (r_flag==1)
	rotatefunc(z, y, nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[i]=z[i];

	for (i=0; i<nx; i++)
	{
		y[i]*=pow(alpha,1.0*i/(nx-1)/2);
	}

	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

    double ret = 0.0;
    for (i=0; i<nx; i++)
    {
        ret += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
    }

	return ret;
}

double step_rastrigin_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Noncontinuous Rastrigin's  */
{
    int i;
	double alpha=10.0,beta=0.2;
	shiftfunc(x, y, nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]=y[i]*5.12/100;
    }

	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

    for (i=0; i<nx; i++)
	{
		if (fabs(z[i])>0.5)
		z[i]=floor(2*z[i]+0.5)/2;
	}

    oszfunc (z, y, nx);
    asyfunc (y, z, nx, beta);

	if (r_flag==1)
	rotatefunc(z, y, nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[i]=z[i];

	for (i=0; i<nx; i++)
	{
		y[i]*=pow(alpha,1.0*i/(nx-1)/2);
	}

	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

    double ret = 0.0;
    for (i=0; i<nx; i++)
    {
        ret += (z[i]*z[i] - 10.0*cos(2.0*PI*z[i]) + 10.0);
    }

	return ret;
}

double schwefel_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Schwefel's  */
{
    int i;
	double tmp;
	shiftfunc(x, y, nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]*=1000/100;
    }
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

	for (i=0; i<nx; i++)
		y[i] = z[i]*pow(10.0,1.0*i/(nx-1)/2.0);

	for (i=0; i<nx; i++)
		z[i] = y[i]+4.209687462275036e+002;
	
    double ret = 0.0;
    for (i=0; i<nx; i++)
	{
		if (z[i]>500)
		{
			ret-=(500.0-fmod(z[i],500))*sin(pow(500.0-fmod(z[i],500),0.5));
			tmp=(z[i]-500.0)/100;
			ret+= tmp*tmp/nx;
		}
		else if (z[i]<-500)
		{
			ret-=(-500.0+fmod(fabs(z[i]),500))*sin(pow(500.0-fmod(fabs(z[i]),500),0.5));
			tmp=(z[i]+500.0)/100;
			ret+= tmp*tmp/nx;
		}
		else
			ret-=z[i]*sin(pow(fabs(z[i]),0.5));
    }
    ret=4.189828872724338e+002*nx+ret;

	return ret;
}

double katsuura_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Katsuura  */
{
    int i,j;
	double temp,tmp1,tmp2,tmp3;
	tmp3=pow(1.0*nx,1.2);
	shiftfunc(x, y, nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]*=5.0/100.0;
    }
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

	for (i=0; i<nx; i++)
		z[i] *=pow(100.0,1.0*i/(nx-1)/2.0);

	if (r_flag==1)
	rotatefunc(z, y, nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[i]=z[i];

    double ret=1.0;
    for (i=0; i<nx; i++)
	{
		temp=0.0;
		for (j=1; j<=32; j++)
		{
			tmp1=pow(2.0,j);
			tmp2=tmp1*y[i];
			temp += fabs(tmp2-floor(tmp2+0.5))/tmp1;
		}
		ret *= pow(1.0+(i+1)*temp,10.0/tmp3);
    }
	tmp1=10.0/nx/nx;
    ret=ret*tmp1-tmp1;

	return ret;
}

double bi_rastrigin_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Lunacek Bi_rastrigin Function */
{
    int i;
	double mu0=2.5,d=1.0,s,mu1,tmp,tmp1,tmp2;
	double *tmpx;
	tmpx=(double *)malloc(sizeof(double)  *  nx);
	s=1.0-1.0/(2.0*pow(nx+20.0,0.5)-8.2);
	mu1=-pow((mu0*mu0-d)/s,0.5);

	shiftfunc(x, y, nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]*=10.0/100.0;
    }

	for (i = 0; i < nx; i++)
    {
		tmpx[i]=2*y[i];
        if (Os[i] < 0.)
            tmpx[i] *= -1.;
    }

	for (i=0; i<nx; i++)
	{
		z[i]=tmpx[i];
		tmpx[i] += mu0;
	}
	if (r_flag==1)
	rotatefunc(z, y, nx, Mr);
	else
    for (i=0; i<nx; i++)
		y[i]=z[i];

	for (i=0; i<nx; i++)
		y[i] *=pow(100.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(y, z, nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

    tmp1=0.0;tmp2=0.0;
    for (i=0; i<nx; i++)
	{
		tmp = tmpx[i]-mu0;
		tmp1 += tmp*tmp;
		tmp = tmpx[i]-mu1;
		tmp2 += tmp*tmp;
    }
	tmp2 *= s;
	tmp2 += d*nx;
	tmp=0;
	for (i=0; i<nx; i++)
	{
		tmp+=cos(2.0*PI*z[i]);
    }
	
	double ret = 0.0;
	if(tmp1<tmp2)
		ret = tmp1;
	else
		ret = tmp2;
	ret += 10.0*(nx-tmp);
	free(tmpx);

	return ret;
}

double grie_rosen_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Griewank-Rosenbrock  */
{
    int i;
    double temp,tmp1,tmp2;

	shiftfunc(x, y, nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[i]=y[i]*5/100;
    }
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

	for (i=0; i<nx; i++)//shift to orgin
    {
        z[i]=y[i]+1;
    }

    double ret=0.0;
    for (i=0; i<nx-1; i++)
    {
		tmp1 = z[i]*z[i]-z[i+1];
		tmp2 = z[i]-1.0;
        temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
        ret += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
	tmp1 = z[nx-1]*z[nx-1]-z[0];
	tmp2 = z[nx-1]-1.0;
    temp = 100.0*tmp1*tmp1 + tmp2*tmp2;;
    ret += (temp*temp)/4000.0 - cos(temp) + 1.0 ;

	return ret;
}

double escaffer6_func (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Expanded Scaffer¡¯s F6  */
{
    int i;
    double temp1, temp2;
	shiftfunc(x, y, nx, Os);
	if (r_flag==1)
	rotatefunc(y, z, nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

	asyfunc (z, y, nx, 0.5);
	if (r_flag==1)
	rotatefunc(y, z, nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		z[i]=y[i];

    double ret = 0.0;
    for (i=0; i<nx-1; i++)
    {
        temp1 = sin(sqrt(z[i]*z[i]+z[i+1]*z[i+1]));
		temp1 =temp1*temp1;
        temp2 = 1.0 + 0.001*(z[i]*z[i]+z[i+1]*z[i+1]);
        ret += 0.5 + (temp1-0.5)/(temp2*temp2);
    }
    temp1 = sin(sqrt(z[nx-1]*z[nx-1]+z[0]*z[0]));
	temp1 =temp1*temp1;
    temp2 = 1.0 + 0.001*(z[nx-1]*z[nx-1]+z[0]*z[0]);
    ret += 0.5 + (temp1-0.5)/(temp2*temp2);

	return ret;
}

double cf01 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 1 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10, 20, 30, 40, 50};
	double bias[5] = {0, 100, 200, 300, 400};
	double ret = 0.0;
	i=0;
	fit[i]=rosenbrock_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+4;
	i=1;
	fit[i]=dif_powers_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+10;
	i=2;
	fit[i]=bent_cigar_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+30;
	i=3;
	fit[i]=discus_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+10;
	i=4;
	fit[i]=sphere_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],0);
	fit[i]=10000*fit[i]/1e+5;
	ret = cf_cal(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf02 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 2 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	for(i=0;i<cf_num;i++)
	{
		fit[i]=schwefel_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	}
	double ret = cf_cal(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf03 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 3 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	for(i=0;i<cf_num;i++)
	{
		fit[i]=schwefel_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	}
	double ret = cf_cal(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf04 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 4 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	i=0;
	fit[i]=schwefel_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	fit[i]=rastrigin_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	fit[i]=weierstrass_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/400;
	double ret = cf_cal(x, nx, Os, delta,bias,fit,cf_num);
	
	return ret;
}

double cf05 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 4 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {10,30,50};
	double bias[3] = {0, 100, 200};
	i=0;
	fit[i]=schwefel_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	fit[i]=rastrigin_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	fit[i]=weierstrass_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/400;
	double ret = cf_cal(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf06 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 6 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,10,10,10,10};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	fit[i]=schwefel_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	fit[i]=rastrigin_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	fit[i]=ellips_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/1e+10;
	i=3;
	fit[i]=weierstrass_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/400;
	i=4;
	fit[i]=griewank_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/100;
	double ret = cf_cal(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf07 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 7 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,10,10,20,20};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	fit[i]=griewank_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/100;
	i=1;
	fit[i]=rastrigin_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+3;
	i=2;
	fit[i]=schwefel_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=3;
	fit[i]=weierstrass_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/400;
	i=4;
	fit[i]=sphere_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],0);
	fit[i]=10000*fit[i]/1e+5;
	double ret = cf_cal(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf08 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 8 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,20,30,40,50};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	fit[i]=grie_rosen_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=1;
	fit[i]=schaffer_F7_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/4e+6;
	i=2;
	fit[i]=schwefel_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=3;
	fit[i]=escaffer6_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/2e+7;
	i=4;
	fit[i]=sphere_func(x,nx,&Os[i*nx],&Mr[i*nx*nx],0);
	fit[i]=10000*fit[i]/1e+5;
	double ret = cf_cal(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cec13::cec13_problems(const double * x, int nx, int func_num) //nx is the dimension
{
	double ret=0;
	//const double * o = new double[dim];
	//const double * M = new double[dim*dim];
	//bool isGood = true;

	switch(func_num)
		{
		case 1:	
			ret=sphere_func(x,nx,OShift,M,0);
			ret+=-1400.0;
			break;
		case 2:	
			ret=ellips_func(x,nx,OShift,M,1);
			ret+=-1300.0;
			break;
		case 3:	
			ret=bent_cigar_func(x,nx,OShift,M,1);
			ret+=-1200.0;
			break;
		case 4:	
			ret=discus_func(x,nx,OShift,M,1);
			ret+=-1100.0;
			break;
		case 5:
			ret=dif_powers_func(x,nx,OShift,M,0);
			ret+=-1000.0;
			break;
		case 6:
			ret=rosenbrock_func(x,nx,OShift,M,1);
			ret+=-900.0;
			break;
		case 7:	
			ret=schaffer_F7_func(x,nx,OShift,M,1);
			ret+=-800.0;
			break;
		case 8:	
			ret=ackley_func(x,nx,OShift,M,1);
			ret+=-700.0;
			break;
		case 9:	
			ret=weierstrass_func(x,nx,OShift,M,1);
			ret+=-600.0;
			break;
		case 10:	
			ret=griewank_func(x,nx,OShift,M,1);
			ret+=-500.0;
			break;
		case 11:	
			ret=rastrigin_func(x,nx,OShift,M,0);
			ret+=-400.0;
			break;
		case 12:	
			ret=rastrigin_func(x,nx,OShift,M,1);
			ret+=-300.0;
			break;
		case 13:	
			ret=step_rastrigin_func(x,nx,OShift,M,1);
			ret+=-200.0;
			break;
		case 14:	
			ret=schwefel_func(x,nx,OShift,M,0);
			ret+=-100.0;
			break;
		case 15:	
			ret=schwefel_func(x,nx,OShift,M,1);
			ret+=100.0;
			break;
		case 16:	
			ret=katsuura_func(x,nx,OShift,M,1);
			ret+=200.0;
			break;
		case 17:	
			ret=bi_rastrigin_func(x,nx,OShift,M,0);
			ret+=300.0;
			break;
		case 18:	
			ret=bi_rastrigin_func(x,nx,OShift,M,1);
			ret+=400.0;
			break;
		case 19:	
			ret=grie_rosen_func(x,nx,OShift,M,1);
			ret+=500.0;
			break;
		case 20:	
			ret=escaffer6_func(x,nx,OShift,M,1);
			ret+=600.0;
			break;
		case 21:	
			ret=cf01(x,nx,OShift,M,1);
			ret+=700.0;
			break;
		case 22:	
			ret=cf02(x,nx,OShift,M,0);
			ret+=800.0;
			break;
		case 23:	
			ret=cf03(x,nx,OShift,M,1);
			ret+=900.0;
			break;
		case 24:	
			ret=cf04(x,nx,OShift,M,1);
			ret+=1000.0;
			break;
		case 25:	
			ret=cf05(x,nx,OShift,M,1);
			ret+=1100.0;
			break;
		case 26:
			ret=cf06(x,nx,OShift,M,1);
			ret+=1200.0;
			break;
		case 27:
			ret=cf07(x,nx,OShift,M,1);
			ret+=1300.0;
			break;
		case 28:
			ret=cf08(x,nx,OShift,M,1);
			ret+=1400.0;
			break;
		}

		/*for (int i=0;i<nx;i++)
		{
			cout<<OShift[i]<<" ";
		}
		getchar();*/
		extern Solution solution;
		for (int i=0;i<nx;i++)
		{
			solution.GlobalOptima[i] = OShift[i];
		}

	return ret;
}