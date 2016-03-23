#include "stdafx.h"
//#include <WINDOWS.H>      
//#include <stdio.h>
//#include <math.h>
//#include <malloc.h>
//#include "cec14.h"

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

using namespace std;

cec14::cec14()
{
	
}

cec14::~cec14()
{
	
}

double *OShift_cec14,*M_cec14,*y_cec14,*z_cec14,*x_bound_cec14;
int ini_flag_cec14,n_flag_cec14,func_flag_cec14,*SS_cec14;


void shiftfunc_cec14 (const double *x, double *xshift, int nx,double *Os)
{
	int i;
    for (i=0; i<nx; i++)
    {
        xshift[i]=x[i]-Os[i];
    }
}

void rotatefunc_cec14 (const double *x, double *xrot, int nx,double *Mr)
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

void sr_func_cec14 (const double *x, double *sr_x, int nx, double *Os,double *Mr, double sh_rate, int s_flag,int r_flag) /* shift and rotate */
{
	int i;
	if (s_flag==1)
	{
		if (r_flag==1)
		{	
			shiftfunc_cec14(x, y_cec14, nx, Os);
			for (i=0; i<nx; i++)//shrink to the original search range
			{
				y_cec14[i]=y_cec14[i]*sh_rate;
			}
			rotatefunc_cec14(y_cec14, sr_x, nx, Mr);
		}
		else
		{
			shiftfunc_cec14(x, sr_x, nx, Os);
			for (i=0; i<nx; i++)//shrink to the original search range
			{
				sr_x[i]=sr_x[i]*sh_rate;
			}
		}
	}
	else
	{	

		if (r_flag==1)
		{	
			for (i=0; i<nx; i++)//shrink to the original search range
			{
				y_cec14[i]=x[i]*sh_rate;
			}
			rotatefunc_cec14(y_cec14, sr_x, nx, Mr);
		}
		else
		for (i=0; i<nx; i++)//shrink to the original search range
		{
			sr_x[i]=x[i]*sh_rate;
		}
	}
}

void asyfunc_cec14 (const double *x, double *xasy, int nx, double beta)
{
	int i;
    for (i=0; i<nx; i++)
    {
		if (x[i]>0)
        xasy[i]=pow(x[i],1.0+beta*i/(nx-1)*pow(x[i],0.5));
    }
}

void oszfunc_cec14 (const double *x, double *xosz, int nx)
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

double cf_cal_cec14(const double *x, int nx, double *Os,double * delta,double * bias,double * fit, int cf_num)
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


double sphere_func_cec14 (const double *x, int nx, double *Os, double *Mr, int s_flag, int r_flag) /* Sphere */
{
	int i;
	double ret = 0.0;
	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */	
	for (i=0; i<nx; i++)
	{					
		ret += z_cec14[i]*z_cec14[i];
	}

	return ret;
}

double ellips_func_cec14 (const double *x, int nx, double *Os,double *Mr, int s_flag, int r_flag) /* Ellipsoidal */
{
    int i;
	double ret = 0.0;
	sr_func_cec14 (x, z_cec14, nx, Os, Mr,1.0, s_flag, r_flag); /* shift and rotate */
	for (i=0; i<nx; i++)
	{
       ret += pow(10.0,6.0*i/(nx-1))*z_cec14[i]*z_cec14[i];
	}

	return ret;
}

double bent_cigar_func_cec14 (const double *x, int nx, double *Os,double *Mr, int s_flag, int r_flag) /* Bent_Cigar */
{
    int i;
	sr_func_cec14 (x, z_cec14, nx, Os, Mr,1.0, s_flag, r_flag); /* shift and rotate */

	double ret = z_cec14[0]*z_cec14[0];
	for (i=1; i<nx; i++)
	{
		ret += pow(10.0,6.0)*z_cec14[i]*z_cec14[i];
	}

	return ret;
}

double discus_func_cec14 (const double *x, int nx, double *Os,double *Mr, int s_flag, int r_flag) /* Discus */
{
    int i;
	sr_func_cec14 (x, z_cec14, nx, Os, Mr,1.0, s_flag, r_flag); /* shift and rotate */
	double ret = pow(10.0,6.0)*z_cec14[0]*z_cec14[0];
	for (i=1; i<nx; i++)
	{
		ret += z_cec14[i]*z_cec14[i];
	}

	return ret;
}

double dif_powers_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Different Powers */
{
	int i;
	double ret = 0.0;
	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		ret += pow(fabs(z_cec14[i]),2+4*i/(nx-1));
	}
	ret=pow(ret,0.5);

	return ret;
}

double rosenbrock_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Rosenbrock's */
{
    int i;
	double tmp1,tmp2;
	double ret = 0.0;
	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 2.048/100.0, s_flag, r_flag); /* shift and rotate */
	z_cec14[0] += 1.0;//shift to orgin
	for (i=0; i<nx-1; i++)
	{
		z_cec14[i+1] += 1.0;//shift to orgin
		tmp1=z_cec14[i]*z_cec14[i]-z_cec14[i+1];
		tmp2=z_cec14[i]-1.0;
		ret += 100.0*tmp1*tmp1 +tmp2*tmp2;
	}

	return ret;
}

double schaffer_F7_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Schwefel's 1.2  */
{
    int i;
	double tmp;
    double ret = 0.0;
	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
	for (i=0; i<nx-1; i++)	
	{
		z_cec14[i]=pow(y_cec14[i]*y_cec14[i]+y_cec14[i+1]*y_cec14[i+1],0.5);
		tmp=sin(50.0*pow(z_cec14[i],0.2));
		ret += pow(z_cec14[i],0.5)+pow(z_cec14[i],0.5)*tmp*tmp ;
	}
	ret = ret*ret/(nx-1)/(nx-1);

	return ret;
}

double ackley_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Ackley's  */
{
    int i;
    double sum1, sum2;
    sum1 = 0.0;
    sum2 = 0.0;

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		sum1 += z_cec14[i]*z_cec14[i];
		sum2 += cos(2.0*PI*z_cec14[i]);
	}
	sum1 = -0.2*sqrt(sum1/nx);
	sum2 /= nx;
	double ret =  E - 20.0*exp(sum1) - exp(sum2) +20.0;

	return ret;
}

double weierstrass_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Weierstrass's  */
{
    int i,j,k_max;
    double sum,sum2, a, b;
    a = 0.5;
    b = 3.0;
    k_max = 20;
    double ret = 0.0;

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 0.5/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		sum = 0.0;
		sum2 = 0.0;
		for (j=0; j<=k_max; j++)
		{
			sum += pow(a,j)*cos(2.0*PI*pow(b,j)*(z_cec14[i]+0.5));
			sum2 += pow(a,j)*cos(2.0*PI*pow(b,j)*0.5);
		}
		ret += sum;
	}
	ret -= nx*sum2;

	return ret;
}

double griewank_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Griewank's  */
{
    int i;
    double s, p;
    s = 0.0;
    p = 1.0;

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 600.0/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		s += z_cec14[i]*z_cec14[i];
		p *= cos(z_cec14[i]/sqrt(1.0+i));
	}
	double ret = 1.0 + s/4000.0 - p;

	return ret;
}

double rastrigin_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Rastrigin's  */
{
    int i;
	double ret = 0.0;

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 5.12/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		ret += (z_cec14[i]*z_cec14[i] - 10.0*cos(2.0*PI*z_cec14[i]) + 10.0);
	}

	return ret;
}

double step_rastrigin_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Noncontinuous Rastrigin's  */
{
    int i;
	double ret=0.0;
	for (i=0; i<nx; i++)
	{
		if (fabs(y_cec14[i]-Os[i])>0.5)
		y_cec14[i]=Os[i]+floor(2*(y_cec14[i]-Os[i])+0.5)/2;
	}

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 5.12/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		ret += (z_cec14[i]*z_cec14[i] - 10.0*cos(2.0*PI*z_cec14[i]) + 10.0);
	}

	return ret;
}

double schwefel_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Schwefel's  */
{
    int i;
	double tmp;
	double ret=0.0;

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1000.0/100.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		z_cec14[i] += 4.209687462275036e+002;
		if (z_cec14[i]>500)
		{
			ret-=(500.0-fmod(z_cec14[i],500))*sin(pow(500.0-fmod(z_cec14[i],500),0.5));
			tmp=(z_cec14[i]-500.0)/100;
			ret+= tmp*tmp/nx;
		}
		else if (z_cec14[i]<-500)
		{
			ret-=(-500.0+fmod(fabs(z_cec14[i]),500))*sin(pow(500.0-fmod(fabs(z_cec14[i]),500),0.5));
			tmp=(z_cec14[i]+500.0)/100;
			ret+= tmp*tmp/nx;
		}
		else
			ret-=z_cec14[i]*sin(pow(fabs(z_cec14[i]),0.5));
	}
	ret +=4.189828872724338e+002*nx;

	return ret;
}

double katsuura_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Katsuura  */
{
    int i,j;
	double temp,tmp1,tmp2,tmp3;
	double ret=1.0;
	tmp3=pow(1.0*nx,1.2);

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */

    for (i=0; i<nx; i++)
	{
		temp=0.0;
		for (j=1; j<=32; j++)
		{
			tmp1=pow(2.0,j);
			tmp2=tmp1*z_cec14[i];
			temp += fabs(tmp2-floor(tmp2+0.5))/tmp1;
		}
		ret *= pow(1.0+(i+1)*temp,10.0/tmp3);
    }
	tmp1=10.0/nx/nx;
    ret=ret*tmp1-tmp1;

	return ret;
}

double bi_rastrigin_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Lunacek Bi_rastrigin Function */
{
    int i;
	double mu0=2.5,d=1.0,s,mu1,tmp,tmp1,tmp2;
	double *tmpx;
	tmpx=(double *)malloc(sizeof(double)  *  nx);
	s=1.0-1.0/(2.0*pow(nx+20.0,0.5)-8.2);
	mu1=-pow((mu0*mu0-d)/s,0.5);
	double ret = 0.0;

	if (s_flag==1)
		shiftfunc_cec14(x, y_cec14, nx, Os);
	else
	{
		for (i=0; i<nx; i++)//shrink to the orginal search range
		{
			y_cec14[i] = x[i];
		}
	}
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y_cec14[i] *= 10.0/100.0;
    }

	for (i = 0; i < nx; i++)
    {
		tmpx[i]=2*y_cec14[i];
        if (Os[i] < 0.0)
            tmpx[i] *= -1.;
    }
	for (i=0; i<nx; i++)
	{
		z_cec14[i]=tmpx[i];
		tmpx[i] += mu0;
	}
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
	tmp=0.0;

	if (r_flag==1)
	{
		rotatefunc_cec14(z_cec14, y_cec14, nx, Mr);
		for (i=0; i<nx; i++)
		{
			tmp+=cos(2.0*PI*y_cec14[i]);
		}	
		if(tmp1<tmp2)
			ret = tmp1;
		else
			ret = tmp2;
		ret += 10.0*(nx-tmp);
	}
	else
	{
		for (i=0; i<nx; i++)
		{
			tmp+=cos(2.0*PI*z_cec14[i]);
		}	
		if(tmp1<tmp2)
			ret = tmp1;
		else
			ret = tmp2;
		ret += 10.0*(nx-tmp);
	}

	free(tmpx);

	return ret;
}

double grie_rosen_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Griewank-Rosenbrock  */
{
    int i;
    double temp,tmp1,tmp2;
    double ret=0.0;

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */

	z_cec14[0] += 1.0;//shift to orgin
    for (i=0; i<nx-1; i++)
    {
		z_cec14[i+1] += 1.0;//shift to orgin
		tmp1 = z_cec14[i]*z_cec14[i]-z_cec14[i+1];
		tmp2 = z_cec14[i]-1.0;
        temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
         ret += (temp*temp)/4000.0 - cos(temp) + 1.0; 
    }
	tmp1 = z_cec14[nx-1]*z_cec14[nx-1]-z_cec14[0];
	tmp2 = z_cec14[nx-1]-1.0;
    temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
    ret += (temp*temp)/4000.0 - cos(temp) + 1.0 ;

	return ret;
}

double escaffer6_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* Expanded Scaffer??s F6  */
{
    int i;
    double temp1, temp2;

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

    double ret = 0.0;
    for (i=0; i<nx-1; i++)
    {
        temp1 = sin(sqrt(z_cec14[i]*z_cec14[i]+z_cec14[i+1]*z_cec14[i+1]));
		temp1 =temp1*temp1;
        temp2 = 1.0 + 0.001*(z_cec14[i]*z_cec14[i]+z_cec14[i+1]*z_cec14[i+1]);
        ret += 0.5 + (temp1-0.5)/(temp2*temp2);
    }
    temp1 = sin(sqrt(z_cec14[nx-1]*z_cec14[nx-1]+z_cec14[0]*z_cec14[0]));
	temp1 =temp1*temp1;
    temp2 = 1.0 + 0.001*(z_cec14[nx-1]*z_cec14[nx-1]+z_cec14[0]*z_cec14[0]);
    ret += 0.5 + (temp1-0.5)/(temp2*temp2);

	return ret;
}

double happycat_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
/* original global optimum: [-1,-1,...,-1] */
{
	int i;
	double alpha,r2,sum_z;
	alpha=1.0/8.0;
	
	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */

	r2 = 0.0;
	sum_z=0.0;
    for (i=0; i<nx; i++)
    {
		z_cec14[i]=z_cec14[i]-1.0;//shift to orgin
        r2 += z_cec14[i]*z_cec14[i];
		sum_z += z_cec14[i];
    }
    double ret=pow(fabs(r2-nx),2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;

	return ret;
}

double hgbat_func_cec14 (const double *x, int nx, double *Os,double *Mr,int s_flag, int r_flag) /* HGBat, provdided by Hans-Georg Beyer (HGB)*/
/* original global optimum: [-1,-1,...,-1] */
{
	int i;
	double alpha,r2,sum_z;
	alpha=1.0/4.0;

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 5.0/100.0, s_flag, r_flag); /* shift and rotate */

	r2 = 0.0;
	sum_z=0.0;
    for (i=0; i<nx; i++)
    {
		z_cec14[i]=z_cec14[i]-1.0;//shift to orgin
        r2 += z_cec14[i]*z_cec14[i];
		sum_z += z_cec14[i];
    }
    double ret=pow(fabs(pow(r2,2.0)-pow(sum_z,2.0)),2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;

	return ret;
}

double hf01_cec14 (const double *x, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 1 */
{
	int i,tmp,cf_num=3;
	double fit[3];
	int G[3],G_nx[3];
	double Gp[3]={0.3,0.3,0.4};

	tmp=0;
	for (i=0; i<cf_num-1; i++)
	{
		G_nx[i] = int(ceil(Gp[i]*nx));
		tmp += G_nx[i];
	}
	G_nx[cf_num-1]=nx-tmp;
	G[0]=0;
	for (i=1; i<cf_num; i++)
	{
		G[i] = G[i-1]+G_nx[i-1];
	}

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		y_cec14[i]=z_cec14[S[i]-1];
	}
	i=0;
	fit[i]=schwefel_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=1;
	fit[i]=rastrigin_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=2;
	fit[i]=ellips_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	double ret=0.0;
	for(i=0;i<cf_num;i++)
	{
		ret += fit[i];
	}

	return ret;
}

double hf02_cec14 (const double *x, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 2 */
{
	int i,tmp,cf_num=3;
	double fit[3];
	int G[3],G_nx[3];
	double Gp[3]={0.3,0.3,0.4};

	tmp=0;
	for (i=0; i<cf_num-1; i++)
	{
		G_nx[i] = int(ceil(Gp[i]*nx));
		tmp += G_nx[i];
	}
	G_nx[cf_num-1]=nx-tmp;

	G[0]=0;
	for (i=1; i<cf_num; i++)
	{
		G[i] = G[i-1]+G_nx[i-1];
	}

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		y_cec14[i]=z_cec14[S[i]-1];
	}
	i=0;
	fit[i]=bent_cigar_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=1;
	fit[i]=hgbat_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=2;
	fit[i]=rastrigin_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);

	double ret=0.0;
	for(i=0;i<cf_num;i++)
	{
		ret += fit[i];
	}

	return ret;
}

double hf03_cec14 (const double *x, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 3 */
{
	int i,tmp,cf_num=4;
	double fit[4];
	int G[4],G_nx[4];
	double Gp[4]={0.2,0.2,0.3,0.3};

	tmp=0;
	for (i=0; i<cf_num-1; i++)
	{
		G_nx[i] = int(ceil(Gp[i]*nx));
		tmp += G_nx[i];
	}
	G_nx[cf_num-1]=nx-tmp;

	G[0]=0;
	for (i=1; i<cf_num; i++)
	{
		G[i] = G[i-1]+G_nx[i-1];
	}

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		y_cec14[i]=z_cec14[S[i]-1];
	}
	i=0;
	fit[i]=griewank_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=1;
	fit[i]=weierstrass_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=2;
	fit[i]=rosenbrock_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=3;
	fit[i]=escaffer6_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	
	double ret=0.0;
	for(i=0;i<cf_num;i++)
	{
		ret += fit[i];
	}

	return ret;
}

double hf04_cec14 (const double *x, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 4 */
{
	int i,tmp,cf_num=4;
	double fit[4];
	int G[4],G_nx[4];
	double Gp[4]={0.2,0.2,0.3,0.3};

	tmp=0;
	for (i=0; i<cf_num-1; i++)
	{
		G_nx[i] = int(ceil(Gp[i]*nx));
		tmp += G_nx[i];
	}
	G_nx[cf_num-1]=nx-tmp;

	G[0]=0;
	for (i=1; i<cf_num; i++)
	{
		G[i] = G[i-1]+G_nx[i-1];
	}

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		y_cec14[i]=z_cec14[S[i]-1];
	}
	i=0;
	fit[i]=hgbat_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=1;
	fit[i]=discus_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=2;
	fit[i]=grie_rosen_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=3;
	fit[i]=rastrigin_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	
	double ret=0.0;
	for(i=0;i<cf_num;i++)
	{
		ret += fit[i];
	}

	return ret;
}

double hf05_cec14 (const double *x, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 5 */
{
	int i,tmp,cf_num=5;
	double fit[5];
	int G[5],G_nx[5];
	double Gp[5]={0.1,0.2,0.2,0.2,0.3};

	tmp=0;
	for (i=0; i<cf_num-1; i++)
	{
		G_nx[i] = int(ceil(Gp[i]*nx));
		tmp += G_nx[i];
	}
	G_nx[cf_num-1]=nx-tmp;

	G[0]=0;
	for (i=1; i<cf_num; i++)
	{
		G[i] = G[i-1]+G_nx[i-1];
	}

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		y_cec14[i]=z_cec14[S[i]-1];
	}
	i=0;
	fit[i]=escaffer6_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0); 
	i=1;
	fit[i]=hgbat_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=2;
	fit[i]=rosenbrock_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=3;
	fit[i]=schwefel_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=4;
	fit[i]=ellips_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);

	double ret=0.0;
	for(i=0;i<cf_num;i++)
	{
		ret += fit[i];
	}

	return ret;
}

double hf06_cec14 (const double *x, int nx, double *Os,double *Mr,int *S,int s_flag,int r_flag) /* Hybrid Function 6 */
{
	int i,tmp,cf_num=5;
	double fit[5];
	int G[5],G_nx[5];
	double Gp[5]={0.1,0.2,0.2,0.2,0.3};

	tmp=0;
	for (i=0; i<cf_num-1; i++)
	{
		G_nx[i] = int(ceil(Gp[i]*nx));
		tmp += G_nx[i];
	}
	G_nx[cf_num-1]=nx-tmp;

	G[0]=0;
	for (i=1; i<cf_num; i++)
	{
		G[i] = G[i-1]+G_nx[i-1];
	}

	sr_func_cec14 (x, z_cec14, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

	for (i=0; i<nx; i++)
	{
		y_cec14[i]=z_cec14[S[i]-1];
	}
	i=0;
	fit[i]=katsuura_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=1;
	fit[i]=happycat_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=2;
	fit[i]=grie_rosen_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=3;
	fit[i]=schwefel_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	i=4;
	fit[i]=ackley_func_cec14(&y_cec14[G[i]],G_nx[i],Os,Mr,0,0);
	
	double ret=0.0;
	for(i=0;i<cf_num;i++)
	{
		ret += fit[i];
	}

	return ret;
}

double cf01_cec14 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 1 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10, 20, 30, 40, 50};
	double bias[5] = {0, 100, 200, 300, 400};
	
	i=0;
	fit[i]=rosenbrock_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1e+4;	
	i=1;
	fit[i]=ellips_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1e+10;
	i=2;
	fit[i]=bent_cigar_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1e+30;	
	i=3;
	fit[i]=discus_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1e+10;	
	i=4;
	fit[i]=ellips_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,0);
	fit[i]=10000*fit[i]/1e+10;	
	
	double ret = cf_cal_cec14(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf02_cec14 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 2 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};

	i=0;
	fit[i]=schwefel_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,0);
	i=1;
	fit[i]=rastrigin_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	i=2;
	fit[i]=hgbat_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	
	double ret = cf_cal_cec14(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf03_cec14 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 3 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {10,30,50};
	double bias[3] = {0, 100, 200};
	i=0;
	fit[i]=schwefel_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	fit[i]=rastrigin_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	fit[i]=ellips_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=1000*fit[i]/1e+10;
	
	double ret = cf_cal_cec14(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf04_cec14 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 4 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,10,10,10,10};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	fit[i]=schwefel_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	fit[i]=happycat_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	fit[i]=ellips_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=1000*fit[i]/1e+10;
	i=3;
	fit[i]=weierstrass_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=1000*fit[i]/400;
	i=4;
	fit[i]=griewank_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=1000*fit[i]/100;
	
	double ret = cf_cal_cec14(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf05_cec14 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 4 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,10,10,20,20};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	fit[i]=hgbat_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1000;
	i=1;
	fit[i]=rastrigin_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1e+3;
	i=2;
	fit[i]=schwefel_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=3;
	fit[i]=weierstrass_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/400;
	i=4;
	fit[i]=ellips_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1e+10;
	
	double ret = cf_cal_cec14(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf06_cec14 (const double *x, int nx, double *Os,double *Mr,int r_flag) /* Composition Function 6 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,20,30,40,50};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	fit[i]=grie_rosen_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=1;
	fit[i]=happycat_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1e+3;
	i=2;
	fit[i]=schwefel_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=3;
	fit[i]=escaffer6_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/2e+7;
	i=4;
	fit[i]=ellips_func_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],1,r_flag);
	fit[i]=10000*fit[i]/1e+10;
	
	double ret = cf_cal_cec14(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf07_cec14 (const double *x, int nx, double *Os,double *Mr,int *SS_cec14,int r_flag) /* Composition Function 7 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {10,30,50};
	double bias[3] = {0, 100, 200};
	i=0;
	fit[i]=hf01_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],&SS_cec14[i*nx],1,r_flag);
	i=1;
	fit[i]=hf02_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],&SS_cec14[i*nx],1,r_flag);
	i=2;
	fit[i]=hf03_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],&SS_cec14[i*nx],1,r_flag);
	
	double ret = cf_cal_cec14(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}

double cf08_cec14 (const double *x, int nx, double *Os,double *Mr,int *SS_cec14,int r_flag) /* Composition Function 8 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {10,30,50};
	double bias[3] = {0, 100, 200};
	i=0;
	fit[i]=hf04_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],&SS_cec14[i*nx],1,r_flag);
	i=1;
	fit[i]=hf05_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],&SS_cec14[i*nx],1,r_flag);
	i=2;
	fit[i]=hf06_cec14(x,nx,&Os[i*nx],&Mr[i*nx*nx],&SS_cec14[i*nx],1,r_flag);
	
	double ret = cf_cal_cec14(x, nx, Os, delta,bias,fit,cf_num);

	return ret;
}


void cec14::FileIO(int dim, int func_num)
{
	int cf_num=10,i,j;
	if (ini_flag_cec14==1)
	{
		if ((n_flag_cec14!=dim)||(func_flag_cec14!=func_num))
		{
			ini_flag_cec14=0;
		}
	}

	if (ini_flag_cec14==0)
	{
		FILE *fpt;
		char FileName[256];
		//free(M_cec14);
		//free(OShift_cec14);
		//free(y_cec14);
		//free(z_cec14);
		//free(x_bound_cec14);
		y_cec14=(double *)malloc(sizeof(double)  *  dim);
		z_cec14=(double *)malloc(sizeof(double)  *  dim);
		x_bound_cec14=(double *)malloc(sizeof(double)  *  dim);
		for (i=0; i<dim; i++)
			x_bound_cec14[i]=100.0;

		if (!(dim==2||dim==10||dim==20||dim==30||dim==50||dim==100))
		{
			printf("\nError: Test functions are only defined for D=2,10,20,30,50,100.\n");
		}
		if (dim==2&&((func_num>=17&&func_num<=22)||(func_num>=29&&func_num<=30)))
		{
			printf("\nError: hf01_cec14,hf02_cec14,hf03_cec14,hf04_cec14,hf05_cec14,hf06_cec14,cf07_cec14&cf08_cec14 are NOT defined for D=2.\n");
		}

		/* Load Matrix M_cec14*/
		sprintf_s(FileName, "input_data_cec14/M_%d_D%d.txt", func_num,dim);
		fopen_s(&fpt,FileName,"r");
		if (fpt==NULL)
		{
		    printf("\n Error: Cannot open input file for reading \n");
		}
		if (func_num<23)
		{
			M_cec14=(double*)malloc(dim*dim*sizeof(double));
			if (M_cec14==NULL)
				printf("\nError: there is insufficient memory available!\n");
			for (i=0; i<dim*dim; i++)
			{
				fscanf_s(fpt,"%Lf",&M_cec14[i]);
			}
		}
		else
		{
			M_cec14=(double*)malloc(cf_num*dim*dim*sizeof(double));
			if (M_cec14==NULL)
				printf("\nError: there is insufficient memory available!\n");
			for (i=0; i<cf_num*dim*dim; i++)
			{
				fscanf_s(fpt,"%Lf",&M_cec14[i]);
			}
		}
		fclose(fpt);
		
		/* Load shift_data */
		sprintf_s(FileName, "input_data_cec14/shift_data_%d.txt", func_num);
		fopen_s(&fpt,FileName,"r");
		if (fpt==NULL)
		{
			printf("\n Error: Cannot open input file for reading \n");
		}

		if (func_num<23)
		{
			OShift_cec14=(double *)malloc(dim*sizeof(double));
			if (OShift_cec14==NULL)
			printf("\nError: there is insufficient memory available!\n");
			for(i=0;i<dim;i++)
			{
				fscanf_s(fpt,"%Lf",&OShift_cec14[i]);
			}
		}
		else
		{
			OShift_cec14=(double *)malloc(dim*cf_num*sizeof(double));
			if (OShift_cec14==NULL)
			printf("\nError: there is insufficient memory available!\n");
			for(i=0;i<cf_num-1;i++)
			{
				for (j=0;j<dim;j++)
				{
					fscanf_s(fpt,"%Lf",&OShift_cec14[i*dim+j]);
				}
				fscanf_s(fpt,"%*[^\n]%*c"); 
			}
			for (j=0;j<dim;j++)
			{
				fscanf_s(fpt,"%Lf",&OShift_cec14[(cf_num-1)*dim+j]);
			}
				
		}
		fclose(fpt);


		/* Load Shuffle_data */
		
		if (func_num>=17&&func_num<=22)
		{
			sprintf_s(FileName, "input_data_cec14/shuffle_data_%d_D%d.txt", func_num, dim);
			fopen_s(&fpt,FileName,"r");
			if (fpt==NULL)
			{
				printf("\n Error: Cannot open input file for reading \n");
			}
			SS_cec14=(int *)malloc(dim*sizeof(int));
			if (SS_cec14==NULL)
				printf("\nError: there is insufficient memory available!\n");
			for(i=0;i<dim;i++)
			{
				fscanf_s(fpt,"%d",&SS_cec14[i]);
			}	
			fclose(fpt);
		}
		else if (func_num==29||func_num==30)
		{
			sprintf_s(FileName, "input_data_cec14/shuffle_data_%d_D%d.txt", func_num, dim);
			fopen_s(&fpt,FileName,"r");
			if (fpt==NULL)
			{
				printf("\n Error: Cannot open input file for reading \n");
			}
			SS_cec14=(int *)malloc(dim*cf_num*sizeof(int));
			if (SS_cec14==NULL)
				printf("\nError: there is insufficient memory available!\n");
			for(i=0;i<dim*cf_num;i++)
			{
				fscanf_s(fpt,"%d",&SS_cec14[i]);
			}
			fclose(fpt);
		}
		

		n_flag_cec14=dim;
		func_flag_cec14=func_num;
		ini_flag_cec14=1;
		//printf("Function has been initialized!\n");
	}
}

void cec14::Release()
{
	free(M_cec14);//Rotation matrix, size n*n
	free(OShift_cec14);//shifting vector, size n
	free(y_cec14);//shifting resultant vector
	//free(z_cec14);
	free(x_bound_cec14);
}

double cec14::cec14_problems(const double * x, int dim, int func_num)
{
	double ret = 0.0;
	
	switch(func_num)
		{
		case 1:	
			ret=ellips_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=100.0;
			break;
		case 2:	
			ret=bent_cigar_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=200.0;
			break;
		case 3:	
			ret=discus_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=300.0;
			break;
		case 4:	
			ret=rosenbrock_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=400.0;
			break;
		case 5:
			ret=ackley_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=500.0;
			break;
		case 6:
			ret=weierstrass_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=600.0;
			break;
		case 7:	
			ret=griewank_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=700.0;
			break;
		case 8:	
			ret=rastrigin_func_cec14(x,dim,OShift_cec14,M_cec14,1,0);
			ret+=800.0;
			break;
		case 9:	
			ret=rastrigin_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=900.0;
			break;
		case 10:	
			ret=schwefel_func_cec14(x,dim,OShift_cec14,M_cec14,1,0);
			ret+=1000.0;
			break;
		case 11:	
			ret=schwefel_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=1100.0;
			break;
		case 12:	
			ret=katsuura_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=1200.0;
			break;
		case 13:	
			ret=happycat_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=1300.0;
			break;
		case 14:	
			ret=hgbat_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=1400.0;
			break;
		case 15:	
			ret=grie_rosen_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=1500.0;
			break;
		case 16:	
			ret=escaffer6_func_cec14(x,dim,OShift_cec14,M_cec14,1,1);
			ret+=1600.0;
			break;
		case 17:	
			ret=hf01_cec14(x,dim,OShift_cec14,M_cec14,SS_cec14,1,1);
			ret+=1700.0;
			break;
		case 18:	
			ret=hf02_cec14(x,dim,OShift_cec14,M_cec14,SS_cec14,1,1);
			ret+=1800.0;
			break;
		case 19:	
			ret=hf03_cec14(x,dim,OShift_cec14,M_cec14,SS_cec14,1,1);
			ret+=1900.0;
			break;
		case 20:	
			ret=hf04_cec14(x,dim,OShift_cec14,M_cec14,SS_cec14,1,1);
			ret+=2000.0;
			break;
		case 21:	
			ret=hf05_cec14(x,dim,OShift_cec14,M_cec14,SS_cec14,1,1);
			ret+=2100.0;
			break;
		case 22:	
			ret=hf06_cec14(x,dim,OShift_cec14,M_cec14,SS_cec14,1,1);
			ret+=2200.0;
			break;
		case 23:	
			ret=cf01_cec14(x,dim,OShift_cec14,M_cec14,1);
			ret+=2300.0;
			break;
		case 24:	
			ret=cf02_cec14(x,dim,OShift_cec14,M_cec14,1);
			ret+=2400.0;
			break;
		case 25:	
			ret=cf03_cec14(x,dim,OShift_cec14,M_cec14,1);
			ret+=2500.0;
			break;
		case 26:
			ret=cf04_cec14(x,dim,OShift_cec14,M_cec14,1);
			ret+=2600.0;
			break;
		case 27:
			ret=cf05_cec14(x,dim,OShift_cec14,M_cec14,1);
			ret+=2700.0;
			break;
		case 28:
			ret=cf06_cec14(x,dim,OShift_cec14,M_cec14,1);
			ret+=2800.0;
			break;
		case 29:
			ret=cf07_cec14(x,dim,OShift_cec14,M_cec14,SS_cec14,1);
			ret+=2900.0;
			break;
		case 30:
			ret=cf08_cec14(x,dim,OShift_cec14,M_cec14,SS_cec14,1);
			ret+=3000.0;
			break;
		default:
			printf("\nError: There are only 30 test functions in this test suite!\n");
			ret = 0.0;
			break;
		}

		extern Solution solution;
		for (int i=0;i<dim;i++)
		{
			solution.GlobalOptima[i] = OShift_cec14[i];
		}

	return ret;
}
