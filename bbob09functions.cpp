#include "stdafx.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <iostream>
//#include "solution.h"
#include "bbob09functions.h"
#include "bbob09supportfunctions.h"

using namespace std;

BBOB09::BBOB09()
{
}

BBOB09::~BBOB09()
{
}

extern Solution solution;
BBOB09support bbob09support;
BBOB_variable bbob09variable;

#define NHIGHPEAKS21 101
#define NHIGHPEAKS22 21

int DIM;
//static unsigned int funcId;
/*double * tmpvect;
double *xshift;
double * tmx;
double ** rotation;
double ** rot2;
double ** linearTF;
double * peaks21;
double * peaks22;
int * rperm21;
int * rperm22;
double ** Xlocal21;
double ** Xlocal22;
double ** arrScales21;
double ** arrScales22;*/

/*
 * Noiseless functions testbed. All functions are ranged in [-5, 5]^DIM.
 */
/*isInitDone status changes when either DIM or trialid change.*/
/*it also changes when a new initialisation has been done (since it rewrites the values of Xopt, Fopt...)*/

void f1_par_cal()
{
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09support.computeXopt(bbob09variable.xshift, rseed);
}
double f1(double* x) {
    /*Sphere function*/
    int i; /*Loop over dim*/
    //static unsigned int funcId = 1;
    double r,res=0;

	//int rseed = solution.Func_num + 10000 * solution.trialID;
    /*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
    //bbob09support.computeXopt(bbob09variable.xshift, rseed);

    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        r = x[i] - bbob09variable.xshift[i];
        res += r * r;
    }
		
	res+=bbob09variable.Fopt;
	//cout<<res<<endl;
    return res;
}

void f2_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift, rseed);
}
double f2(double* x) {
    /* separable ellipsoid with monotone transformation, condition 1e6*/
    int i; /*Loop over dim*/
    //static unsigned int funcId = 2;
    double res=0;
    static double condition = 1e6;
	
	//rseed = solution.Func_num + 10000 * solution.trialID;
	
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	//bbob09support.computeXopt(bbob09variable.xshift, rseed);

    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = x[i] - bbob09variable.xshift[i];
    }
	
    bbob09support.monotoneTFosc(bbob09variable.tmx);
	
    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        res += pow(condition, ((double)i)/((double)(solution.D-1))) * bbob09variable.tmx[i] * bbob09variable.tmx[i];
    }
	
	res+=bbob09variable.Fopt;
    return res;
}

void f3_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift, rseed);
}
double f3(double* x) {
    /* Rastrigin with monotone transformation separable "condition" 10*/
    int i; /*Loop over dim*/
    static unsigned int funcId = 3;
    static double condition = 10.;
    static double beta = 0.2;
    double tmp, tmp2;
    double res=0;
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed);

    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = x[i] - bbob09variable.xshift[i];
    }

    bbob09support.monotoneTFosc(bbob09variable.tmx);
    for (i = 0; i < solution.D; i++)
    {
        tmp = ((double)i)/((double)(solution.D-1));
        if (bbob09variable.tmx[i] > 0)
            bbob09variable.tmx[i] = pow(bbob09variable.tmx[i], 1 + beta * tmp * sqrt(bbob09variable.tmx[i]));
        bbob09variable.tmx[i] = pow(sqrt(condition), tmp) * bbob09variable.tmx[i];
    }
    /* COMPUTATION core*/
    tmp = 0.;
    tmp2 = 0.;
    for (i = 0; i < solution.D; i++)
    {
        tmp += cos(2*M_PI*bbob09variable.tmx[i]);
        tmp2 += bbob09variable.tmx[i]*bbob09variable.tmx[i];
    }
	res = 10 * (solution.D - tmp) + tmp2;
	res+=bbob09variable.Fopt;
    
    return res;
}

void f4_par_cal()
{
	int rseed = 3 + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift, rseed);
	for (int i = 0; i < solution.D; i += 2)
		bbob09variable.xshift[i] = fabs(bbob09variable.xshift[i]); /*Skew*/
}
double f4(double* x) {
    /* skew Rastrigin-Bueche, condition 10, skew-"condition" 100*/
    int i; /*Loop over dim*/
    static unsigned int funcId = 4;
    static double condition = 10.;
    static double alpha = 100.;
    double tmp, tmp2;
    double res=0;

	//rseed = funcId + 10000 * solution.trialID; /* Not the same as before.*/
	/*INITIALIZATION*/
	double Fpen=0.;
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed);
	//for (i = 0; i < solution.D; i += 2)
	//	bbob09variable.xshift[i] = fabs(bbob09variable.xshift[i]); /*Skew*/

    for (i = 0; i < solution.D; i++) {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
            Fpen += tmp * tmp;
    }
    Fpen *= 1e2;
	double tmFopt = bbob09variable.Fopt;
    tmFopt += Fpen;

    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = x[i] - bbob09variable.xshift[i];
    }

    bbob09support.monotoneTFosc(bbob09variable.tmx);
    for (i = 0; i < solution.D; i++)
    {
        if (i % 2 == 0 && bbob09variable.tmx[i] > 0)
            bbob09variable.tmx[i] = sqrt(alpha) * bbob09variable.tmx[i];
        bbob09variable.tmx[i] = pow(sqrt(condition), ((double)i)/((double)(solution.D-1))) * bbob09variable.tmx[i];
    }
    /* COMPUTATION core*/
    tmp = 0.;
    tmp2 = 0.;
    for (i = 0; i < solution.D; i++)
    {
        tmp += cos(2*M_PI*bbob09variable.tmx[i]);
        tmp2 += bbob09variable.tmx[i]*bbob09variable.tmx[i];
    }
	res = 10 * (solution.D - tmp) + tmp2;
    res+=tmFopt;
    
    return res;
}

void f5_par_cal()
{
	static double alpha = 100.;
	double tmp=0;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift, rseed);
	for (int i = 0; i < solution.D; i ++)
    {
        tmp = pow(sqrt(alpha), ((double)i)/((double)(solution.D-1)));
        if (bbob09variable.xshift[i] > 0)
        {
            bbob09variable.xshift[i] = 5.;
        }
        else if (bbob09variable.xshift[i] < 0)
        {
            bbob09variable.xshift[i] = -5.;
        }
		bbob09variable.Fopt += 5. * tmp;
    }
}
double f5(double* x) {
    /* linear slope*/
    int i;//, rseed; /*Loop over dim*/
    //static unsigned int funcId = 5;
    static double alpha = 100.;
    static double Fadd; /*Treatment is different from other functions.*/
    double Ftrue = 0.;
    double res=0;

    //rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
    //bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift, rseed);
    /*for (i = 0; i < solution.D; i ++)
    {
        tmp = pow(sqrt(alpha), ((double)i)/((double)(solution.D-1)));
        if (bbob09variable.xshift[i] > 0)
        {
            bbob09variable.xshift[i] = 5.;
        }
        else if (bbob09variable.xshift[i] < 0)
        {
            bbob09variable.xshift[i] = -5.;
        }
		bbob09variable.Fopt += 5. * tmp;
    }*/

    /* BOUNDARY HANDLING*/
    /* move "too" good coordinates back into domain*/
    for (i = 0; i < solution.D; i++) {
        if ((bbob09variable.xshift[i] == 5.) && (x[i] > 5))
            bbob09variable.tmx[i] = 5.;
        else if ((bbob09variable.xshift[i] == -5.) && (x[i] < -5))
            bbob09variable.tmx[i] = -5.;
        else
            bbob09variable.tmx[i] = x[i];
    }

    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        if (bbob09variable.xshift[i] > 0) {
            res -= pow(sqrt(alpha), ((double)i)/((double)(solution.D-1))) * bbob09variable.tmx[i];
        } else {
            res += pow(sqrt(alpha), ((double)i)/((double)(solution.D-1))) * bbob09variable.tmx[i];
        }
    }
    res += bbob09variable.Fopt;

    return res;
}

void f6_par_cal()
{
	static double condition = 10.;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift, rseed);

	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
    bbob09support.computeRotation(bbob09variable.rot2, rseed);
    /* decouple scaling from function definition*/
    for (int i = 0; i < solution.D; i ++)
    {
        for (int j = 0; j < solution.D; j++)
        {
            bbob09variable.linearTF[i][j] = 0.;
            for (int k = 0; k < solution.D; k++) 
			{
                bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(solution.D-1))) * bbob09variable.rot2[k][j];
            }
        }
    }
}
double f6(double* x) {
    /* attractive sector function*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 6;
    static double alpha = 100.;
    double Ftrue = 0.;
    static double condition = 10.;
    //rseed = funcId + 10000 * solution.trialID;
    /*INITIALIZATION*/
    //bbob09variable.Fopt = bbob09support.computeFopt(funcId);
    //bbob09support.computeXopt(bbob09variable.xshift,rseed);
	
    //bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
    //bbob09support.computeRotation(bbob09variable.rot2, rseed);
    /* decouple scaling from function definition*/
    /*for (i = 0; i < solution.D; i ++)
    {
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.linearTF[i][j] = 0.;
            for (k = 0; k < solution.D; k++) 
			{
                bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(solution.D-1))) * bbob09variable.rot2[k][j];
            }
        }
    }*/

    /* BOUNDARY HANDLING*/
    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
	{
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmx[i] += bbob09variable.linearTF[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        if (bbob09variable.tmx[i] * bbob09variable.xshift[i] > 0)
            bbob09variable.tmx[i] *= alpha;
        Ftrue += bbob09variable.tmx[i] * bbob09variable.tmx[i];
    }
	
    /*MonotoneTFosc...*/
    if (Ftrue > 0)
    {
        Ftrue = pow(exp(log(Ftrue)/0.1 + 0.49*(sin(log(Ftrue)/0.1) + sin(0.79*log(Ftrue)/0.1))), 0.1);
    }
    else if (Ftrue < 0)
    {
        Ftrue = -pow(exp(log(-Ftrue)/0.1 + 0.49*(sin(0.55 * log(-Ftrue)/0.1) + sin(0.31*log(-Ftrue)/0.1))), 0.1);
    }
    Ftrue = pow(Ftrue, 0.9);
    Ftrue += bbob09variable.Fopt;

    return Ftrue;
}

void f7_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift, rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
    bbob09support.computeRotation(bbob09variable.rot2, rseed);    
}
double f7(double* x) {
    /* step-ellipsoid, condition 100*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 7;
    static double condition = 100.;
    static double alpha = 10.;
    double x1, tmp, Fpen = 0., Ftrue = 0.;

	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed);
	//bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	//bbob09support.computeRotation(bbob09variable.rot2, rseed);

    /* BOUNDARY HANDLING*/
    for (i = 0; i < solution.D; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
    tmFopt += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++) 
	{
        bbob09variable.tmpvect[i] = 0.;
        tmp = sqrt(pow(condition/10., ((double)i)/((double)(solution.D-1))));
        for (j = 0; j < solution.D; j++)
		{
            bbob09variable.tmpvect[i] += tmp * bbob09variable.rot2[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
    }
    x1 = bbob09variable.tmpvect[0];

    for (i = 0; i < solution.D; i++) {
        if (fabs(bbob09variable.tmpvect[i]) > 0.5)
            bbob09variable.tmpvect[i] = bbob09support.round(bbob09variable.tmpvect[i]);
        else
            bbob09variable.tmpvect[i] = bbob09support.round(alpha * bbob09variable.tmpvect[i])/alpha;
    }

    for (i = 0; i < solution.D; i++) {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmx[i] += bbob09variable.rotation[i][j] * bbob09variable.tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(solution.D-1))) * bbob09variable.tmx[i] * bbob09variable.tmx[i];
    }
    Ftrue = 0.1 * bbob09support.fmax(1e-4 * fabs(x1), Ftrue);

    Ftrue += tmFopt;

    return Ftrue;
}

void f8_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift, rseed);
	for (int i = 0; i < solution.D; i ++)
        bbob09variable.xshift[i] *= 0.75;  
}
double f8(double* x) {
    /* Rosenbrock, non-rotated*/
    static unsigned int funcId = 8;
    int i; /*Loop over dim*/
    double tmp, Ftrue = 0.;
	static double scales = bbob09support.fmax(1., sqrt(solution.D) / 8.);
	scales = bbob09support.fmax(1., sqrt(solution.D) / 8.);
	//rseed = funcId + 10000 * solution.trialID;
	//scales = bbob09support.fmax(1., sqrt(solution.D) / 8.);
	/*INITIALIZATION*/
    //bbob09variable.Fopt = bbob09support.computeFopt(funcId);
    //bbob09support.computeXopt(bbob09variable.xshift,rseed);
    //for (i = 0; i < solution.D; i ++)
     //   bbob09variable.xshift[i] *= 0.75;

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++) {
        bbob09variable.tmx[i] = scales * (x[i] - bbob09variable.xshift[i]) + 1;
    }

    /* COMPUTATION core*/
    for (i = 0; i < solution.D - 1; i++)
    {
        tmp = (bbob09variable.tmx[i] * bbob09variable.tmx[i] - bbob09variable.tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < solution.D - 1; i ++)
    {
        tmp = (bbob09variable.tmx[i] - 1.);
        Ftrue += tmp * tmp;
    }
    Ftrue += bbob09variable.Fopt;

    return Ftrue;
}

void f9_par_cal()
{
	double scales = 0;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeRotation(bbob09variable.rotation, rseed);
	scales = bbob09support.fmax(1., sqrt(solution.D) / 8.);
	for (int i = 0; i < solution.D; i ++)
	{
        for (int j = 0; j < solution.D; j++)
            bbob09variable.linearTF[i][j] = scales * bbob09variable.rotation[i][j];
    }
}
double f9(double* x)
{
	/* Rosenbrock, rotated*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 9;
    double tmp, Ftrue = 0.;

	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	/* computeXopt(rseed, DIM);*/
	/*bbob09support.computeRotation(bbob09variable.rotation, rseed);
	scales = bbob09support.fmax(1., sqrt(solution.D) / 8.);
	for (i = 0; i < solution.D; i ++)
	{
        for (j = 0; j < solution.D; j++)
            bbob09variable.linearTF[i][j] = scales * bbob09variable.rotation[i][j];
    }*/

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++) {
        bbob09variable.tmx[i] = 0.5;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmx[i] += bbob09variable.linearTF[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < solution.D - 1; i++)
    {
        tmp = (bbob09variable.tmx[i] * bbob09variable.tmx[i] - bbob09variable.tmx[i+1]);
        Ftrue += tmp * tmp;
    }
    Ftrue *= 1e2;
    for (i = 0; i < solution.D - 1; i ++)
    {
       tmp = (bbob09variable.tmx[i] - 1.);
        Ftrue += tmp * tmp;
    }

    Ftrue += bbob09variable.Fopt;
    
    return Ftrue;
}

void f10_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
}
double f10(double* x)
{
	/* ellipsoid with monotone transformation, condition 1e6*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 10;
    static double condition = 1e6;
    double Ftrue = 0.;
    
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed);
	//bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmx[i] += bbob09variable.rotation[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
    }

    bbob09support.monotoneTFosc(bbob09variable.tmx);
    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        Ftrue += pow(condition, ((double)i)/((double)(solution.D-1))) * bbob09variable.tmx[i] * bbob09variable.tmx[i];
    }
    Ftrue += bbob09variable.Fopt;

    return Ftrue;
}

void f11_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
}
double f11(double* x)
{
	/* discus (tablet) with monotone transformation, condition 1e6*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 11;
    static double condition = 1e6;
    double Ftrue;
    
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed);
	//bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmx[i] += bbob09variable.rotation[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
    }

    bbob09support.monotoneTFosc(bbob09variable.tmx);

    /* COMPUTATION core*/
    Ftrue = condition * bbob09variable.tmx[0] * bbob09variable.tmx[0];
    for (i = 1; i < solution.D; i++)
    {
        Ftrue += bbob09variable.tmx[i] * bbob09variable.tmx[i];
    }
    Ftrue += bbob09variable.Fopt;
    
    return Ftrue;
}

void f12_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
}
double f12(double* x)
{
	/* bent cigar with asymmetric space distortion, condition 1e6*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 12;
    static double condition = 1e6;
    static double beta = 0.5;
    double Ftrue;

	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed + 1000000);
	//bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmpvect[i] = 0.;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmpvect[i] += bbob09variable.rotation[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
        if (bbob09variable.tmpvect[i] > 0)
        {
            bbob09variable.tmpvect[i] = pow(bbob09variable.tmpvect[i], 1 + beta * ((double)i)/((double)(solution.D-1)) * sqrt(bbob09variable.tmpvect[i]));
        }
    }

    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmx[i] += bbob09variable.rotation[i][j] * bbob09variable.tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    Ftrue = bbob09variable.tmx[0] * bbob09variable.tmx[0];
    for (i = 1; i < solution.D; i++)
    {
        Ftrue += condition * bbob09variable.tmx[i] * bbob09variable.tmx[i];
    }
    Ftrue += bbob09variable.Fopt;
    
    return Ftrue;
}

void f13_par_cal()
{
	static double condition = 10.;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);

	for (int i=0;i<solution.D;i++)
	{
		for (int j = 0; j < solution.D; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (int k = 0; k < solution.D; k++)
			{
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(solution.D-1))) * bbob09variable.rot2[k][j];
			}
		}
	}
}
double f13(double* x)
{
	/* sharp ridge*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 13;
    static double condition = 10.;
    static double alpha = 100.;
    double Ftrue = 0.;

	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	/*bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);

	for (i=0;i<solution.D;i++)
	{
		for (j = 0; j < solution.D; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (k = 0; k < solution.D; k++)
			{
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(solution.D-1))) * bbob09variable.rot2[k][j];
			}
		}
	}*/

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmx[i] += bbob09variable.linearTF[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 1; i < solution.D; i++)
    {
        Ftrue += bbob09variable.tmx[i] * bbob09variable.tmx[i];
    }
    Ftrue = alpha * sqrt(Ftrue);
    Ftrue += bbob09variable.tmx[0] * bbob09variable.tmx[0];

    Ftrue += bbob09variable.Fopt;
    
    return Ftrue;
}

void f14_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
}
double f14(double* x)
{
	/* sum of different powers, between x^2 and x^6*/
    int i, j;
    static unsigned int funcId = 14;
    static double alpha = 4.;
    double Ftrue = 0.;
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed);
	//bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++) {
            bbob09variable.tmx[i] += bbob09variable.rotation[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        Ftrue += pow(fabs(bbob09variable.tmx[i]), 2. + alpha * ((double)i)/((double)(solution.D-1)));
    }
    Ftrue = sqrt(Ftrue);

    Ftrue += bbob09variable.Fopt;

    return Ftrue;
}

void f15_par_cal()
{
	static double condition = 10.;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
	for (int i = 0; i < solution.D; i++)
	{
		for (int j = 0; j < solution.D; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (int k = 0; k < solution.D; k++) {
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(solution.D-1))) * bbob09variable.rot2[k][j];
			}
		}
	}
}
double f15(double* x)
{
	/* Rastrigin with asymmetric non-linear distortion, "condition" 10*/
	int i,j;
    static unsigned int funcId = 15;
    static double condition = 10.;
    static double beta = 0.2;
    double tmp = 0., tmp2 = 0., Ftrue;

	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	/*bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
	for (i = 0; i < solution.D; i++)
	{
		for (j = 0; j < solution.D; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (k = 0; k < solution.D; k++) {
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(solution.D-1))) * bbob09variable.rot2[k][j];
			}
		}
	}*/

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmpvect[i] = 0.;
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.tmpvect[i] += bbob09variable.rotation[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
    }

    bbob09support.monotoneTFosc(bbob09variable.tmpvect);
    for (i = 0; i < solution.D; i++)
    {
        if (bbob09variable.tmpvect[i] > 0)
            bbob09variable.tmpvect[i] = pow(bbob09variable.tmpvect[i], 1 + beta * ((double)i)/((double)(solution.D-1)) * sqrt(bbob09variable.tmpvect[i]));
    }
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.tmx[i] += bbob09variable.linearTF[i][j] * bbob09variable.tmpvect[j];
        }
    }
    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        tmp += cos(2. * M_PI * bbob09variable.tmx[i]);
        tmp2 += bbob09variable.tmx[i] * bbob09variable.tmx[i];
    }
    Ftrue = 10. * ((double)solution.D - tmp) + tmp2;
    Ftrue += bbob09variable.Fopt;
    
    return Ftrue;
}

void f16_par_cal()
{
	static double condition = 100.;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
	for (int i = 0; i < solution.D; i++)
	{
		for (int j = 0; j < solution.D; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (int k = 0; k < solution.D; k++) {
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(1./sqrt(condition), ((double)k)/((double)(solution.D-1))) * bbob09variable.rot2[k][j];
			}
		}
	}
	bbob09variable.F0 = 0.;
	for (int i = 0; i < 12; i ++) /* number of summands, 20 in CEC2005, 10/12 saves 30% of time*/
	{
		bbob09variable.aK[i] = pow(0.5, (double)i);
		bbob09variable.bK[i] = pow(3., (double)i);
		bbob09variable.F0 += bbob09variable.aK[i] * cos(2 * M_PI * bbob09variable.bK[i] * 0.5);
	}
}
double f16(double* x)
{
	/* Weierstrass, condition 100*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 16;
    static double condition = 100.;
    //static double aK[12];
    //static double bK[12];
    //static double F0;
    double tmp, Fpen = 0., Ftrue = 0.;
    
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	/*bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
	for (i = 0; i < solution.D; i++)
	{
		for (j = 0; j < solution.D; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (k = 0; k < solution.D; k++) {
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(1./sqrt(condition), ((double)k)/((double)(solution.D-1))) * bbob09variable.rot2[k][j];
			}
		}
	}
	bbob09variable.F0 = 0.;*/
	//for (i = 0; i < 12; i ++) /* number of summands, 20 in CEC2005, 10/12 saves 30% of time*/
	/*{
		bbob09variable.aK[i] = pow(0.5, (double)i);
		bbob09variable.bK[i] = pow(3., (double)i);
		bbob09variable.F0 += bbob09variable.aK[i] * cos(2 * M_PI * bbob09variable.bK[i] * 0.5);
	}*/

    /* BOUNDARY HANDLING*/
    for (i = 0; i < solution.D; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
    tmFopt += 10./(double)solution.D * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmpvect[i] = 0.;
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.tmpvect[i] += bbob09variable.rotation[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
    }

    bbob09support.monotoneTFosc(bbob09variable.tmpvect);
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.tmx[i] += bbob09variable.linearTF[i][j] * bbob09variable.tmpvect[j];
        }
    }
    /* COMPUTATION core*/
    for (i = 0; i < solution.D; i++)
    {
        tmp = 0.;
        for (j = 0; j < 12; j++)
        {
            tmp += cos(2 * M_PI * (bbob09variable.tmx[i] + 0.5) * bbob09variable.bK[j]) * bbob09variable.aK[j];
        }
        Ftrue += tmp;
    }
    Ftrue = 10. * pow(Ftrue/(double)solution.D - bbob09variable.F0, 3.);
    Ftrue += tmFopt;
    
    return Ftrue;
}

void f17_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
}
double f17(double* x)
{
	/* Schaffers F7 with asymmetric non-linear transformation, condition 10*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 17;
    static double condition = 10.;
    static double beta = 0.5;
    double tmp, Fpen = 0., Ftrue = 0.;
	
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed);
	//bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	//bbob09support.computeRotation(bbob09variable.rot2, rseed);

    /* BOUNDARY HANDLING*/
    for (i = 0; i < solution.D; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
    tmFopt += 10. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmpvect[i] = 0.;
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.tmpvect[i] += bbob09variable.rotation[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
        if (bbob09variable.tmpvect[i] > 0)
            bbob09variable.tmpvect[i] = pow(bbob09variable.tmpvect[i], 1 + beta * ((double)i)/((double)(solution.D-1)) * sqrt(bbob09variable.tmpvect[i]));
    }

    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        tmp = pow(sqrt(condition), ((double)i)/((double)(solution.D-1)));
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.tmx[i] += tmp * bbob09variable.rot2[i][j] * bbob09variable.tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < solution.D - 1; i++)
    {
        tmp = bbob09variable.tmx[i] * bbob09variable.tmx[i] + bbob09variable.tmx[i+1] * bbob09variable.tmx[i+1];
        Ftrue += pow(tmp, 0.25) * (pow(sin(50 * pow(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = pow(Ftrue/(double)(solution.D - 1), 2.);
    Ftrue += tmFopt;
    
    return Ftrue;
}

void f18_par_cal()
{
	int rseed = 17 + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
}
double f18(double* x)
{
	/* Schaffers F7 with asymmetric non-linear transformation, condition 1000*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 18;
    static double condition = 1e3;
    static double beta = 0.5;
    double tmp, Fpen = 0., Ftrue = 0.;
	
	//rseed = 17 + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.computeXopt(bbob09variable.xshift,rseed);
	//bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	//bbob09support.computeRotation(bbob09variable.rot2, rseed);

    /* BOUNDARY HANDLING*/
    for (i = 0; i < solution.D; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
    tmFopt += 10. * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmpvect[i] = 0.;
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.tmpvect[i] += bbob09variable.rotation[i][j] * (x[j] - bbob09variable.xshift[j]);
        }
        if (bbob09variable.tmpvect[i] > 0)
            bbob09variable.tmpvect[i] = pow(bbob09variable.tmpvect[i], 1. + beta * ((double)i)/((double)(solution.D-1)) * sqrt(bbob09variable.tmpvect[i]));
    }

    for (i = 0; i < solution.D; i++)
    {
        bbob09variable.tmx[i] = 0.;
        tmp = pow(sqrt(condition), ((double)i)/((double)(solution.D-1)));
        for (j = 0; j < solution.D; j++)
        {
            bbob09variable.tmx[i] += tmp * bbob09variable.rot2[i][j] * bbob09variable.tmpvect[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < solution.D - 1; i++)
    {
        tmp = bbob09variable.tmx[i] * bbob09variable.tmx[i] + bbob09variable.tmx[i+1] * bbob09variable.tmx[i+1];
        Ftrue += pow(tmp, 0.25) * (pow(sin(50. * pow(tmp, 0.1)), 2.) + 1.);
    }
    Ftrue = pow(Ftrue/(double)(solution.D - 1), 2.);
    Ftrue += tmFopt;
    
    return Ftrue;
}

void f19_par_cal()
{
	double scales = 0;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	scales = bbob09support.fmax(1., sqrt(solution.D) / 8.);
	bbob09support.computeRotation(bbob09variable.rotation, rseed);
	for (int i = 0; i < DIM; i ++)
	{
		for (int j = 0; j < DIM; j++)
		{
			bbob09variable.linearTF[i][j] = scales * bbob09variable.rotation[i][j];
		}
	}
	for (int i = 0; i < DIM; i++)
	{
		bbob09variable.xshift[i] = 0.;
		cout<<bbob09variable.xshift[i]<<" ";
		for (int j = 0; j < DIM; j++)
		{
			bbob09variable.xshift[i] += bbob09variable.linearTF[j][i] * 0.5/scales/scales;
		}
	}
}
double f19(double* x)
{
	/* F8F2 sum of Griewank-Rosenbrock 2-D blocks*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 19;
    double F2, tmp = 0., tmp2, Ftrue = 0.;
	
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	/*bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	scales = bbob09support.fmax(1., sqrt(solution.D) / 8.);
	bbob09support.computeRotation(bbob09variable.rotation, rseed);
	for (i = 0; i < DIM; i ++)
	{
		for (j = 0; j < DIM; j++)
		{
			bbob09variable.linearTF[i][j] = scales * bbob09variable.rotation[i][j];
		}
	}
	for (i = 0; i < DIM; i++)
	{
		bbob09variable.xshift[i] = 0.;
		for (j = 0; j < DIM; j++)
		{
			bbob09variable.xshift[i] += bbob09variable.linearTF[j][i] * 0.5/scales/scales;
		}
	}*/

    /* BOUNDARY HANDLING*/

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++) {
        bbob09variable.tmx[i] = 0.5;
        for (j = 0; j < DIM; j++) {
            bbob09variable.tmx[i] += bbob09variable.linearTF[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < DIM - 1; i++)
    {
        tmp2 = bbob09variable.tmx[i] * bbob09variable.tmx[i] -bbob09variable.tmx[i+1];
        F2 = 100. * tmp2 * tmp2;
        tmp2 = 1 - bbob09variable.tmx[i];
        F2 += tmp2 * tmp2;
        tmp += F2 / 4000. - cos(F2);
    }
    Ftrue = 10. + 10. * tmp / (double)(DIM - 1);

    Ftrue += bbob09variable.Fopt;
    
    return Ftrue;
}

void f20_par_cal()
{
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.unif(bbob09variable.tmpvect, DIM, rseed);
	for (int i = 0; i < DIM; i++)
	{
		bbob09variable.xshift[i] = 0.5 * 4.2096874633;
		if (bbob09variable.tmpvect[i] - 0.5 < 0)
			bbob09variable.xshift[i] *= -1.;
	}
}
double f20(double* x)
{
	/* Schwefel with tridiagonal variable transformation*/
    int i; /*Loop over dim*/
    static unsigned int funcId = 20;
    static double condition = 10.;
    double tmp, Fpen = 0., Ftrue = 0.;
	
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	//bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	//bbob09support.unif(bbob09variable.tmpvect, DIM, rseed);
	//for (i = 0; i < DIM; i++)
	//{
	//	bbob09variable.xshift[i] = 0.5 * 4.2096874633;
	//	if (bbob09variable.tmpvect[i] - 0.5 < 0)
	//		bbob09variable.xshift[i] *= -1.;
	//}

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        bbob09variable.tmpvect[i] = 2. * x[i];
        if (bbob09variable.xshift[i] < 0.)
            bbob09variable.tmpvect[i] *= -1.;
    }

    bbob09variable.tmx[0] = bbob09variable.tmpvect[0];
    for (i = 1; i < DIM; i++)
    {
        bbob09variable.tmx[i] = bbob09variable.tmpvect[i] + 0.25 * (bbob09variable.tmpvect[i-1] - 2. * fabs(bbob09variable.xshift[i-1]));
    }

    for (i = 0; i < DIM; i++)
    {
        bbob09variable.tmx[i] -= 2 * fabs(bbob09variable.xshift[i]);
        bbob09variable.tmx[i] *= pow(sqrt(condition), ((double)i)/((double)(DIM-1)));
        bbob09variable.tmx[i] = 100. * (bbob09variable.tmx[i] + 2 * fabs(bbob09variable.xshift[i]));
    }

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabs(bbob09variable.tmx[i]) - 500.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
    tmFopt += 0.01 * Fpen;

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        Ftrue += bbob09variable.tmx[i] * sin(sqrt(fabs(bbob09variable.tmx[i])));
    }
    Ftrue = 0.01 * ((418.9828872724339) - Ftrue / (double)DIM);

	Ftrue += tmFopt;
    
    return Ftrue;
}

int compare_doubles21(const void *a, const void *b)
{
    double temp = bbob09variable.peaks21[*(const int*)a] - bbob09variable.peaks21[*(const int*)b];
    if (temp > 0)
        return 1;
    else if (temp < 0)
        return -1;
    else
        return 0;
}
void f21_par_cal()
{
	static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS21];
	static double fitvalues[2] = {1.1, 9.1};
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeRotation(bbob09variable.rotation, rseed);
	bbob09support.unif(bbob09variable.peaks21, NHIGHPEAKS21 - 1, rseed);
	for (int i = 0; i < NHIGHPEAKS21 - 1; i++)
		bbob09variable.rperm21[i] = i;
	qsort(bbob09variable.rperm21, NHIGHPEAKS21 - 1, sizeof(int), compare_doubles21);

	/* Random permutation*/
	arrCondition[0] = sqrt(maxcondition);
	bbob09variable.peakvalues21[0] = 10;
	for (int i = 1; i < NHIGHPEAKS21; i++)
	{
		arrCondition[i] = pow(maxcondition, (double)(bbob09variable.rperm21[i-1])/((double)(NHIGHPEAKS21-2)));
		bbob09variable.peakvalues21[i] = (double)(i-1)/(double)(NHIGHPEAKS21-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
	}
	
	for (int i = 0; i < NHIGHPEAKS21; i++)
	{
		bbob09support.unif(bbob09variable.peaks21, DIM, rseed + 1000 * i);
		for (int j = 0; j < DIM; j++)
			bbob09variable.rperm21[j] = j;
		qsort(bbob09variable.rperm21, DIM, sizeof(int), compare_doubles21);
		for (int j = 0; j < DIM; j++)
		{
			bbob09variable.arrScales21[i][j] = pow(arrCondition[i], ((double)bbob09variable.rperm21[j])/((double)(DIM-1)) - 0.5);
		}
	}

	bbob09support.unif(bbob09variable.peaks21, DIM * NHIGHPEAKS21, rseed);
	
	for (int i = 0; i < DIM; i++)
	{
		bbob09variable.xshift[i] = 0.8 * (10. * bbob09variable.peaks21[i] -5.);
		for (int j = 0; j < NHIGHPEAKS21; j++)
		{
			bbob09variable.Xlocal21[i][j] = 0.;
			for (int k = 0; k < DIM; k++)
			{
				bbob09variable.Xlocal21[i][j] += bbob09variable.rotation[i][k] * (10. * bbob09variable.peaks21[j * DIM + k] -5.);
			}
			if (j == 0)
				bbob09variable.Xlocal21[i][j] *= 0.8;
		}
	}
}
double f21(double* x)
{
	/* Gallagher with 101 Gaussian peaks, condition up to 1000, one global rotation*/
    int i, j; /*Loop over dim*/
    //static unsigned int funcId = 21;
    //static double fitvalues[2] = {1.1, 9.1};
    //static double maxcondition = 1000.;
    //static double arrCondition[NHIGHPEAKS21];
    //static double peakvalues[NHIGHPEAKS21];
    static double a = 0.1;
    double tmp2, f = 0., tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)DIM;
	
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	/*bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	bbob09support.computeRotation(bbob09variable.rotation, rseed);
	bbob09support.unif(bbob09variable.peaks21, NHIGHPEAKS21 - 1, rseed);
	for (i = 0; i < NHIGHPEAKS21 - 1; i++)
		bbob09variable.rperm21[i] = i;
	qsort(bbob09variable.rperm21, NHIGHPEAKS21 - 1, sizeof(int), compare_doubles21);*/

	/* Random permutation*/
	/*arrCondition[0] = sqrt(maxcondition);
	peakvalues[0] = 10;
	for (i = 1; i < NHIGHPEAKS21; i++)
	{
		arrCondition[i] = pow(maxcondition, (double)(bbob09variable.rperm21[i-1])/((double)(NHIGHPEAKS21-2)));
		peakvalues[i] = (double)(i-1)/(double)(NHIGHPEAKS21-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
	}
	
	for (i = 0; i < NHIGHPEAKS21; i++)
	{
		bbob09support.unif(bbob09variable.peaks21, DIM, rseed + 1000 * i);
		for (j = 0; j < DIM; j++)
			bbob09variable.rperm21[j] = j;
		qsort(bbob09variable.rperm21, DIM, sizeof(int), compare_doubles21);
		for (j = 0; j < DIM; j++)
		{
			bbob09variable.arrScales21[i][j] = pow(arrCondition[i], ((double)bbob09variable.rperm21[j])/((double)(DIM-1)) - 0.5);
		}
	}

	bbob09support.unif(bbob09variable.peaks21, DIM * NHIGHPEAKS21, rseed);
	
	for (i = 0; i < DIM; i++)
	{
		bbob09variable.xshift[i] = 0.8 * (10. * bbob09variable.peaks21[i] -5.);
		for (j = 0; j < NHIGHPEAKS21; j++)
		{
			bbob09variable.Xlocal21[i][j] = 0.;
			for (k = 0; k < DIM; k++)
			{
				bbob09variable.Xlocal21[i][j] += bbob09variable.rotation[i][k] * (10. * bbob09variable.peaks21[j * DIM + k] -5.);
			}
			if (j == 0)
				bbob09variable.Xlocal21[i][j] *= 0.8;
		}
	}*/

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
	tmFopt += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            bbob09variable.tmx[i] += bbob09variable.rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS21; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmp = (bbob09variable.tmx[j] - bbob09variable.Xlocal21[j][i]);
            tmp2 += bbob09variable.arrScales21[i][j] * tmp * tmp;
        }
        tmp2 = bbob09variable.peakvalues21[i] * exp(fac * tmp2);
        f = bbob09support.fmax(f, tmp2);
    }

    f = 10. - f;
    if (f > 0)
    {
        Ftrue = log(f)/a;
        Ftrue = pow(exp(Ftrue + 0.49*(sin(Ftrue) + sin(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = log(-f)/a;
        Ftrue = -pow(exp(Ftrue + 0.49*(sin(0.55 * Ftrue) + sin(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;
	Ftrue += tmFopt;
    
    return Ftrue;
}

int compare_doubles22(const void *a, const void *b)
{
    double temp = bbob09variable.peaks22[*(const int*)a] - bbob09variable.peaks22[*(const int*)b];
    if (temp > 0)
        return 1;
    else if (temp < 0)
        return -1;
    else
        return 0;
}
void f22_par_cal()
{
	static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS22];
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeRotation(bbob09variable.rotation, rseed);
	bbob09support.unif(bbob09variable.peaks22, NHIGHPEAKS22 - 1, rseed);
	for (int i = 0; i < NHIGHPEAKS22 - 1; i++)
		bbob09variable.rperm22[i] = i;
	qsort(bbob09variable.rperm22, NHIGHPEAKS22 - 1, sizeof(int), compare_doubles22);
	arrCondition[0] = maxcondition;
	bbob09variable.peakvalues22[0] = 10;
	for (int i = 1; i < NHIGHPEAKS22; i++)
	{
		arrCondition[i] = pow(maxcondition, (double)(bbob09variable.rperm22[i-1])/((double)(NHIGHPEAKS22-2)));
		bbob09variable.peakvalues22[i] = (double)(i-1)/(double)(NHIGHPEAKS22-2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
	}
	for (int i = 0; i < NHIGHPEAKS22; i++)
	{
		bbob09support.unif(bbob09variable.peaks22, DIM, rseed + 1000 * i);
		for (int j = 0; j < DIM; j++)
			bbob09variable.rperm22[j] = j;
		qsort(bbob09variable.rperm22, DIM, sizeof(int), compare_doubles22);
		for (int j = 0; j < DIM; j++)
		{
			bbob09variable.arrScales22[i][j] = pow(arrCondition[i], ((double)bbob09variable.rperm22[j])/((double)(DIM-1)) - 0.5);
		}
	}

	bbob09support.unif(bbob09variable.peaks22, DIM * NHIGHPEAKS22, rseed);
	for (int i = 0; i < DIM; i++)
	{
		bbob09variable.xshift[i] = 0.8 * (9.8 * bbob09variable.peaks22[i] -4.9);
		for (int j = 0; j < NHIGHPEAKS22; j++)
		{
			bbob09variable.Xlocal22[i][j] = 0.;
			for (int k = 0; k < DIM; k++)
			{
				bbob09variable.Xlocal22[i][j] += bbob09variable.rotation[i][k] * (9.8 * bbob09variable.peaks22[j * DIM + k] -4.9);
			}
			if (j == 0)
				bbob09variable.Xlocal22[i][j] *= 0.8;
		}
	}
}
double f22(double* x)
{
	/* Gallagher with 21 Gaussian peaks, condition up to 1000, one global rotation*/
    int i, j; /*Loop over dim*/
    /*static unsigned int funcId = 22;
    static double fitvalues[2] = {1.1, 9.1};
    static double maxcondition = 1000.;
    static double arrCondition[NHIGHPEAKS22];
    static double peakvalues[NHIGHPEAKS22];*/
    static double a = 0.1;
    double tmp2, f = 0., tmp, Fpen = 0., Ftrue = 0.;
    double fac = -0.5 / (double)DIM;
    
    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
    tmFopt += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        bbob09variable.tmx[i] = 0.;
        for (j = 0; j < DIM; j++)
        {
            bbob09variable.tmx[i] += bbob09variable.rotation[i][j] * x[j];
        }
    }

    /* COMPUTATION core*/
    for (i = 0; i < NHIGHPEAKS22; i++)
    {
        tmp2 = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmp = (bbob09variable.tmx[j] - bbob09variable.Xlocal22[j][i]);
            tmp2 += bbob09variable.arrScales22[i][j] * tmp * tmp;
        }
        tmp2 = bbob09variable.peakvalues22[i] * exp(fac * tmp2);
        f = bbob09support.fmax(f, tmp2);
    }

    f = 10. - f;
    if (f > 0)
    {
        Ftrue = log(f)/a;
        Ftrue = pow(exp(Ftrue + 0.49*(sin(Ftrue) + sin(0.79*Ftrue))), a);
    }
    else if (f < 0)
    {
        Ftrue = log(-f)/a;
        Ftrue = -pow(exp(Ftrue + 0.49*(sin(0.55 * Ftrue) + sin(0.31*Ftrue))), a);
    }
    else
        Ftrue = f;

    Ftrue *= Ftrue;
    Ftrue += tmFopt;
    
    return Ftrue;
}

void f23_par_cal()
{
	static double condition = 100.;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
	
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (int k = 0; k < DIM; k++)
			{
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/(double)(DIM - 1)) * bbob09variable.rot2[k][j];
			}
		}
	}
}
double f23(double* x)
{
	/* Katsuura function*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 23;
    static double condition = 100.;
    double Fpen = 0., tmp, Ftrue = 0., arr, prod = 1., tmp2;
    double *ptmx, *plinTF, *ptmp;
	
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	/*bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	bbob09support.computeXopt(bbob09variable.xshift,rseed);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
	
	for (i = 0; i < DIM; i++)
	{
		for (j = 0; j < DIM; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (k = 0; k < DIM; k++)
			{
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/(double)(DIM - 1)) * bbob09variable.rot2[k][j];
			}
		}
	}*/

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
	tmFopt += Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    /* write rotated difference vector into tmx*/
    for (j = 0; j < DIM; j++)  /* store difference vector*/
        bbob09variable.tmpvect[j] = x[j] - bbob09variable.xshift[j];
    for (i = 0; i < DIM; i++) 
	{
        bbob09variable.tmx[i] = 0.;
        ptmx = &bbob09variable.tmx[i];
        plinTF = bbob09variable.linearTF[i];
        ptmp = bbob09variable.tmpvect;
        for (j = 0; j < DIM; j++) 
		{
            *ptmx += *plinTF++ * *ptmp++;
        }
    }

/*     for (i = 0; i < DIM; i++) {
           tmx[i] = 0.;
           for (j = 0; j < DIM; j++) {
               tmx[i] += linearTF[i][j] * (x[j] - Xopt[j]);
           }
       }*/

    /* COMPUTATION core*/
    for (i = 0; i < DIM; i++)
    {
        tmp = 0.;
        for (j = 1; j < 33; j++)
        {
            tmp2 = pow(2., (double)j);
            arr = bbob09variable.tmx[i] * tmp2;
            tmp += fabs(arr - bbob09support.round(arr)) / tmp2;
        }
        tmp = 1. + tmp * (double)(i + 1);
        prod *= tmp;
    }
    Ftrue = 10./(double)DIM/(double)DIM * (-1. + pow(prod, 10./pow((double)DIM, 1.2)));
	Ftrue += tmFopt;
    
	//delete []ptmp;
	//delete []ptmx;
	//delete []plinTF;
    return Ftrue;
}

void f24_par_cal()
{
	static double condition = 100.;
	static double mu1 = 2.5;
	int rseed = solution.Func_num + 10000 * solution.trialID;
	bbob09variable.Fopt = bbob09support.computeFopt(solution.Func_num);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
	bbob09support.gauss(bbob09variable.tmpvect, DIM, rseed);
	for (int i = 0; i < DIM; i++)
	{
		bbob09variable.xshift[i] = 0.5 * mu1;
		if (bbob09variable.tmpvect[i] < 0.)
			bbob09variable.xshift[i] *= -1.;
	}
	
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (int k = 0; k < DIM; k++) 
			{
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(DIM-1))) * bbob09variable.rot2[k][j];
			}
		}
	}
}
double f24(double* x)
{
	/* Lunacek bi-Rastrigin, condition 100*/
    /* in PPSN 2008, Rastrigin part rotated and scaled*/
    int i, j; /*Loop over dim*/
    static unsigned int funcId = 24;
    static double condition = 100.;
    static double mu1 = 2.5;
    double Fpen = 0., tmp, Ftrue = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
    double s = 1. - 0.5 / (sqrt((double)(DIM + 20)) - 4.1);
    static double d = 1.;
    double mu2 = -sqrt((mu1 * mu1 - d) / s);
	
	//rseed = funcId + 10000 * solution.trialID;
	/*INITIALIZATION*/
	/*bbob09variable.Fopt = bbob09support.computeFopt(funcId);
	bbob09support.computeRotation(bbob09variable.rotation, rseed + 1000000);
	bbob09support.computeRotation(bbob09variable.rot2, rseed);
	bbob09support.gauss(bbob09variable.tmpvect, DIM, rseed);
	for (i = 0; i < DIM; i++)
	{
		bbob09variable.xshift[i] = 0.5 * mu1;
		if (bbob09variable.tmpvect[i] < 0.)
			bbob09variable.xshift[i] *= -1.;
	}
	
	for (i = 0; i < DIM; i++)
	{
		for (j = 0; j < DIM; j++)
		{
			bbob09variable.linearTF[i][j] = 0.;
			for (k = 0; k < DIM; k++) 
			{
				bbob09variable.linearTF[i][j] += bbob09variable.rotation[i][k] * pow(sqrt(condition), ((double)k)/((double)(DIM-1))) * bbob09variable.rot2[k][j];
			}
		}
	}*/

    /* BOUNDARY HANDLING*/
    for (i = 0; i < DIM; i++)
    {
        tmp = fabs(x[i]) - 5.;
        if (tmp > 0.)
        {
            Fpen += tmp * tmp;
        }
    }
	double tmFopt = bbob09variable.Fopt;
	tmFopt += 1e4 * Fpen;

    /* TRANSFORMATION IN SEARCH SPACE*/
    for (i = 0; i < DIM; i++)
    {
        bbob09variable.tmx[i] = 2. * x[i];
        if (bbob09variable.xshift[i] < 0.)
            bbob09variable.tmx[i] *= -1.;
    }

    /* COMPUTATION core*/
    tmp = 0.;
    for (i = 0; i < DIM; i++)
    {
        tmp2 += (bbob09variable.tmx[i] - mu1) * (bbob09variable.tmx[i] - mu1);
        tmp3 += (bbob09variable.tmx[i] - mu2) * (bbob09variable.tmx[i] - mu2);
        tmp4 = 0.;
        for (j = 0; j < DIM; j++)
        {
            tmp4 += bbob09variable.linearTF[i][j] * (bbob09variable.tmx[j] - mu1);
        }
        tmp += cos(2 * M_PI * tmp4);
    }
    Ftrue = bbob09support.fmin(tmp2, d * (double)DIM + s * tmp3) + 10. * ((double)DIM - tmp);
	Ftrue += tmFopt;
    
    return Ftrue;
}

void BBOB09::BBOBparameterInitial()
{
	bbob09variable.Fopt = 0;
	DIM = solution.D;
	//funcId=solution.Func_num;
	bbob09variable.xshift = new double[solution.D];
	bbob09variable.rotation = new double*[solution.D];
	bbob09variable.rot2 = new double*[solution.D];
	bbob09variable.linearTF = new double*[solution.D];
	bbob09variable.Xlocal21 = new double*[solution.D];
	bbob09variable.Xlocal22 = new double*[solution.D];
	for (int i=0;i<solution.D;i++)
	{
		bbob09variable.rotation[i]=new double[solution.D];
		bbob09variable.rot2[i]=new double[solution.D];
		bbob09variable.linearTF[i]=new double[solution.D];
		bbob09variable.Xlocal21[i]=new double[NHIGHPEAKS21];
		bbob09variable.Xlocal22[i]=new double[NHIGHPEAKS22];
	}
	bbob09variable.tmx = new double[solution.D];
	bbob09variable.tmpvect = new double[solution.D];
	bbob09variable.peaks21 = new double[solution.D*NHIGHPEAKS21];
	bbob09variable.peaks22 = new double[solution.D*NHIGHPEAKS22];
	int nr21 = int(bbob09support.fmax(DIM, NHIGHPEAKS21 - 1));
	int nr22 = int(bbob09support.fmax(DIM, NHIGHPEAKS22 - 1));
	bbob09variable.rperm21 = new int[nr21];
	bbob09variable.rperm22 = new int[nr22];
	bbob09variable.arrScales21 = new double*[NHIGHPEAKS21];
	bbob09variable.arrScales22 = new double*[NHIGHPEAKS22];
	for (int i=0;i<NHIGHPEAKS21;i++)
	{
		bbob09variable.arrScales21[i]=new double[solution.D];
	}
	for (int i=0;i<NHIGHPEAKS22;i++)
	{
		bbob09variable.arrScales22[i]=new double[solution.D];
	}
	bbob09variable.F0=0;
	for (int i=0;i<12;i++)
	{
		bbob09variable.aK[i]=0;
		bbob09variable.bK[i]=0;
	}
	bbob09variable.peakvalues21 = new double[NHIGHPEAKS21];
	bbob09variable.peakvalues22 = new double[NHIGHPEAKS22];
}

void BBOB09::Release()
{
	delete []bbob09variable.xshift;
	for (int i=0;i<solution.D;i++)
	{
		delete []bbob09variable.rotation[i];
		delete []bbob09variable.rot2[i];
		delete []bbob09variable.linearTF[i];
		delete []bbob09variable.Xlocal21[i];
		delete []bbob09variable.Xlocal22[i];
	}
	delete []bbob09variable.rotation;
	delete []bbob09variable.rot2;
	delete []bbob09variable.linearTF;
	delete []bbob09variable.Xlocal21;
	delete []bbob09variable.Xlocal22;
	delete []bbob09variable.tmx;
	delete []bbob09variable.tmpvect;
	delete []bbob09variable.rperm21;
	delete []bbob09variable.rperm22;
	delete []bbob09variable.peaks21;
	delete []bbob09variable.peaks22;
	for (int i=0;i<NHIGHPEAKS21;i++)
	{
		delete []bbob09variable.arrScales21[i];
	}
	for (int i=0;i<NHIGHPEAKS22;i++)
	{
		delete []bbob09variable.arrScales22[i];
	}
	delete []bbob09variable.arrScales21;
	delete []bbob09variable.arrScales22;
	delete []bbob09variable.peakvalues21;
	delete []bbob09variable.peakvalues22;
}

void BBOB09::CalculateParameter(int no)
{
	switch (no) {
	case 1:
		f1_par_cal();
		break;
	case 2:
		f2_par_cal();
		break;
	case 3:
		f3_par_cal();
		break;
	case 4:
		f4_par_cal();
		break;
	case 5:
		f5_par_cal();
		break;
	case 6:
		f6_par_cal();
		break;
	case 7:
		f7_par_cal();
		break;
	case 8:
		f8_par_cal();
		break;
	case 9:
		f9_par_cal();
		break;
	case 10:
		f10_par_cal();
		break;
	case 11:
		f11_par_cal();
		break;
	case 12:
		f12_par_cal();
		break;
	case 13:
		f13_par_cal();
		break;
	case 14:
		f14_par_cal();
		break;
	case 15:
		f15_par_cal();
		break;
	case 16:
		f16_par_cal();
		break;
	case 17:
		f17_par_cal();
		break;
	case 18:
		f18_par_cal();
		break;
	case 19:
		f19_par_cal();
		break;
	case 20:
		f20_par_cal();
		break;
	case 21:
		f21_par_cal();
		break;
	case 22:
		f22_par_cal();
		break;
	case 23:
		f23_par_cal();
		break;
	case 24:
		f24_par_cal();
		break;
	}
}

double BBOB09::FunctionCalculation(double * x, int no) 
{
	double ret=0;

	if (no < 1 || no > 24) {
		//isGood = false;
		cout<<"Wrong Function Number! Should be 1-24";
		return 0.0;
	}

	switch (no) {
	case 1:
		ret=f1(x);
		break;
	case 2:
		ret=f2(x);
		break;
	case 3:
		ret=f3(x);
		break;
	case 4:
		ret=f4(x);
		break;
	case 5:
		ret=f5(x);
		break;
	case 6:
		ret=f6(x);
		break;
	case 7:
		ret=f7(x);
		break;
	case 8:
		ret=f8(x);
		break;
	case 9:
		ret=f9(x);
		break;
	case 10:
		ret=f10(x);
		break;
	case 11:
		ret=f11(x);
		break;
	case 12:
		ret=f12(x);
		break;
	case 13:
		ret=f13(x);
		break;
	case 14:
		ret=f14(x);
		break;
	case 15:
		ret=f15(x);
		break;
	case 16:
		ret=f16(x);
		break;
	case 17:
		ret=f17(x);
		break;
	case 18:
		ret=f18(x);
		break;
	case 19:
		ret=f19(x);
		break;
	case 20:
		ret=f20(x);
		break;
	case 21:
		ret=f21(x);
		break;
	case 22:
		ret=f22(x);
		break;
	case 23:
		ret=f23(x);
		break;
	case 24:
		ret=f24(x);
		break;
	}
	
	return ret;
}


