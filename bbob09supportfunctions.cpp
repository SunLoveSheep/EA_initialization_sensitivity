#include "stdafx.h"
#include <stdarg.h>
#include <string.h>
#include "solution.h"
#include "bbob09supportfunctions.h"

#define TOL 1e-8
extern Solution solution;
static int seed = -1;
static int seedn = -1;

extern double * peaks;

BBOB09support::BBOB09support()
{
}

BBOB09support::~BBOB09support()
{
}

void unif(double* r, int N, int inseed)
{
    /* generates N uniform numbers with starting seed*/
    int aktseed;
    int tmp;
    int rgrand[32];
    int aktrand;
    int i;

    if (inseed < 0)
        inseed = -inseed;
    if (inseed < 1)
        inseed = 1;
    aktseed = inseed;

    for (i = 39; i >= 0; i--)
    {
        tmp = floor((double)aktseed/(double)127773);
        aktseed = 16807  * (aktseed - tmp * 127773) - 2836 * tmp;
        if (aktseed < 0)
            aktseed = aktseed + 2147483647;
        if (i < 32)
           rgrand[i] = aktseed;
    }
    aktrand = rgrand[0];

    for (i = 0; i < N; i++)
    {
        tmp = floor((double)aktseed/(double)127773);
        aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
        if (aktseed < 0)
            aktseed = aktseed + 2147483647;
        tmp = floor((double)aktrand / (double)67108865);
        aktrand = rgrand[tmp];
        rgrand[tmp] = aktseed;
        r[i] = (double)aktrand/2.147483647e9;
        if (r[i] == 0.)
        {
            printf("Warning: zero sampled(?), set to 1e-99.\n");
            r[i] = 1e-99;
        }
    }
    return;
}

double fmax(double a, double b) {
 return b > a ? b : a;
}

/* set the seed for the noise.
 * If the seeds are larger than 1e9 they are set back to 1 in randn and myrand.
 */
void setNoiseSeed(unsigned int _seed, unsigned int _seedn)
{
    seed = _seed;
    seedn = _seedn;
}

void gauss(double * g, int N, int seed)
{
    /* samples N standard normally distributed numbers
       being the same for a given seed.*/
    int i;
	double *uniftmp = (double*)malloc(sizeof(double) * 2 * solution.D * solution.D);

    unif(uniftmp, 2*N, seed);

    for (i = 0; i < N; i++)
    {
        g[i] = sqrt(-2*log(uniftmp[i])) * cos(2*M_PI*uniftmp[N+i]);
        if (g[i] == 0.)
            g[i] = 1e-99;
    }
	free(uniftmp);
    return;
}

void reshape(double** B, double* vector, int m, int n)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            B[i][j] = vector[j * m + i];
        }
    }
}

double myrand() {
    /*Adaptation of myrand*/
    if (seed == -1)
        seed = time(NULL) % 1000000000; /* cannot be larger than 1e9 */
	double *uniftmp = (double*)malloc(sizeof(double) * 2 * solution.D * solution.D);
    seed ++;
    if (seed > 1e9)
        seed = 1;
    unif(uniftmp, 1, seed);
	double result=uniftmp[0];
	free(uniftmp);
    return result;
}

double randn() {
    /*Adaptation of myrandn*/
    if (seedn == -1)
        seedn = time(NULL) % 1000000000; /* cannot be larger than 1e9 */
	double *uniftmp = (double*)malloc(sizeof(double) * 2 * solution.D * solution.D);
    seedn ++;
    if (seedn > 1e9)
        seedn = 1;
    gauss(uniftmp, 1, seedn);
	double result=uniftmp[0];
	free(uniftmp);
    return result;
}

double FGauss(double Ftrue, double beta)
{
    double Fval = Ftrue * exp(beta * randn());
    Fval += 1.01 * TOL;
    if (Ftrue < TOL) {
        Fval = Ftrue;
    }
    return Fval;
}

double FUniform(double Ftrue, double alpha, double beta)
{
    double Fval = pow(myrand(), beta) * Ftrue * fmax(1., pow(1e9/(Ftrue+1e-99), alpha * myrand()));
    Fval += 1.01 * TOL;
    if (Ftrue < TOL) {
        Fval = Ftrue;
    }
    return Fval;
}

double FCauchy(double Ftrue, double alpha, double p)
{
    double Fval;
    double tmp = randn()/fabs(randn()+1e-199);
    /* tmp is so as to actually do the calls to randn in order for the number
     * of calls to be the same as in the Matlab code.
     */
    if (myrand() < p)
        Fval = Ftrue + alpha * fmax(0., 1e3 + tmp);
    else
        Fval = Ftrue + alpha * 1e3;

    Fval += 1.01 * TOL;
    if (Ftrue < TOL) {
        Fval = Ftrue;
    }
    return Fval;
}

double BBOB09support::round(double a) {
 return floor(a + 0.5);
}
double BBOB09support::fmin(double a, double b) {
 return b < a ? b : a;
}
double BBOB09support::fmax(double a, double b) {
 return b > a ? b : a;
}

void BBOB09support::gauss(double * g, int N, int seed)
{
    /* samples N standard normally distributed numbers
       being the same for a given seed.*/
    int i;
	double *uniftmp = (double*)malloc(sizeof(double) * 2 * solution.D * solution.D);

    unif(uniftmp, 2*N, seed);

    for (i = 0; i < N; i++)
    {
        g[i] = sqrt(-2*log(uniftmp[i])) * cos(2*M_PI*uniftmp[N+i]);
        if (g[i] == 0.)
            g[i] = 1e-99;
    }
	free(uniftmp);
    return;
}

void BBOB09support::unif(double* r, int N, int inseed)
{
    /* generates N uniform numbers with starting seed*/
    int aktseed;
    int tmp;
    int rgrand[32];
    int aktrand;
    int i;

    if (inseed < 0)
        inseed = -inseed;
    if (inseed < 1)
        inseed = 1;
    aktseed = inseed;

    for (i = 39; i >= 0; i--)
    {
        tmp = floor((double)aktseed/(double)127773);
        aktseed = 16807  * (aktseed - tmp * 127773) - 2836 * tmp;
        if (aktseed < 0)
            aktseed = aktseed + 2147483647;
        if (i < 32)
           rgrand[i] = aktseed;
    }
    aktrand = rgrand[0];

    for (i = 0; i < N; i++)
    {
        tmp = floor((double)aktseed/(double)127773);
        aktseed = 16807 * (aktseed - tmp * 127773) - 2836 * tmp;
        if (aktseed < 0)
            aktseed = aktseed + 2147483647;
        tmp = floor((double)aktrand / (double)67108865);
        aktrand = rgrand[tmp];
        rgrand[tmp] = aktseed;
        r[i] = (double)aktrand/2.147483647e9;
        if (r[i] == 0.)
        {
            printf("Warning: zero sampled(?), set to 1e-99.\n");
            r[i] = 1e-99;
        }
    }
    return;
}

void BBOB09support::computeXopt(double *xshift, int seed) {
    int i;
	double *tmpvect = (double*)malloc(sizeof(double) * solution.D);
    unif(tmpvect, solution.D, seed);
    for (i = 0; i < solution.D; i++)
    {
        xshift[i] = 8 * floor(1e4 * tmpvect[i])/1e4 - 4;
        if (xshift[i] == 0.0)
            xshift[i] = -1e-5;
    }
	free(tmpvect);
}

double BBOB09support::computeFopt(int _funcId) {
    int rseed, rrseed;
	double *gval = (double*)malloc(sizeof(double) * 1);
    double *gval2 = (double*)malloc(sizeof(double) * 1);

    if (_funcId == 4)
        rseed = 3;
    else if (_funcId == 18)
        rseed = 17;
    else
        rseed = _funcId;

	rrseed = rseed + 10000 * solution.trialID;
    gauss(gval, 1, rrseed);
    gauss(gval2, 1, rrseed + 1);
	double g1=gval[0];
	double g2=gval2[0];
	free(gval);
	free(gval2);
    return fmin(1000., fmax(-1000., (round(100.*100.*g1/g2)/100.)));
}

void BBOB09support::monotoneTFosc(double* f) {
    double a = 0.1;
    int i;
    for (i = 0; i < solution.D; i++)
    {
        if (f[i] > 0)
        {
            f[i] = log(f[i])/a;
            f[i] = pow(exp(f[i] + 0.49*(sin(f[i]) + sin(0.79*f[i]))), a);
        }
        else if (f[i] < 0)
        {
            f[i] = log(-f[i])/a;
            f[i] = -pow(exp(f[i] + 0.49*(sin(0.55 * f[i]) + sin(0.31*f[i]))), a);
        }
    }
}

void BBOB09support::computeRotation(double ** B, int seed)
{
    double prod;
    int i, j, k; /*Loop over pairs of column vectors*/
	double *gvect = (double*)malloc(sizeof(double) * solution.D * solution.D);
    gauss(gvect, solution.D * solution.D, seed);
    reshape(B, gvect, solution.D, solution.D);
    /*1st coordinate is row, 2nd is column.*/

    for (i = 0; i < solution.D; i++)
    {
        for (j = 0; j < i; j++)
        {
            prod = 0;
            for (k = 0; k < solution.D; k++)
            {
                prod += B[k][i] * B[k][j];
            }
            for (k = 0; k < solution.D; k++)
            {
                B[k][i] -= prod * B[k][j];
            }
        }
        prod = 0;
        for (k = 0; k < solution.D; k++)
        {
            prod += B[k][i] * B[k][i];
        }
        for (k = 0; k < solution.D; k++)
        {
            B[k][i] /= sqrt(prod);
        }
    }
	free(gvect);
}
