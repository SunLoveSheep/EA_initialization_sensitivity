#include "stdafx.h"
#include "DE.h"
#include "ConstraintHandle.h"
/*#include "solution.h"
#include "Initialization.h"
#include "bbob09functions.h"
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "ConstraintHandle.h"
#include "FEoutput.h"
#include <iostream>
#include <string>*/

using std::string;
using namespace std;

DEparameter DEp;

DE::DE()
{
	DEp.F=0.8;
	DEp.CR=0.9;
}

DE::~DE()
{
	
}

//use v=x1+F(x2-x3) to find new parent
void DE::DEmutation(double *X, int *R)
{
	extern Solution solution;
	ConstraintHandle constrainthandle;
	
	int *tempX=new int[solution.EANP];
	
	int r,r1,r2,r3;//xr0=xr1+F(xr2-xr3)
	int rcount=0;
	//to generate four none-repeated random number from 0 to population size
	for (int i=0;i<solution.EANP;i++)
	{
		tempX[i]=i;
	}

	int swapnumber=0,temp=0;//for swapping the two numbers
	for (int i=solution.EANP-1;i>solution.EANP-5;i--)
	{
		swapnumber=rand()%(i+1);
		temp=tempX[i];
		tempX[i]=tempX[swapnumber];
		tempX[swapnumber]=temp;
	}
	
	r=tempX[solution.EANP-1];
	r1=tempX[solution.EANP-2];
	r2=tempX[solution.EANP-3];
	r3=tempX[solution.EANP-4];
	
	//DE mutation
	for (int i=0;i<solution.D;i++)
	{
		X[i]=solution.S[r1][i]+DEp.F*(solution.S[r2][i]-solution.S[r3][i]);
	}
	constrainthandle.Handleconstaint(X);

	R[0]=r;
	delete []tempX;
}

void DE::DEcrossover(double *X, int *R)
{
	extern Solution solution;
	Initialization initialization;
	
	int j=rand()%(solution.D);
	double rn;//a random number generated bewtween 0,1 to compare with CR
	for (int i=0;i<solution.D;i++)
	{
		rn=initialization.RandomGen(0,1);
		if ((rn<=DEp.CR)||(i==j))
		{
			X[i]=X[i];
		}
		else if((rn>DEp.CR)&&(i!=j))
		{
			X[i]=solution.S[R[0]][i];
		}
	}
}

void DE::DEselection(double *X, int *R)
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	double *temp=new double[solution.D];
	for (int i=0;i<solution.D;i++)
	{
		temp[i]=solution.S[R[0]][i];
	}
	double FVorigin=0, FVresult=0;//to store the function value of two compared solution candidates
	//FVorigin=CEC14.cec14_eotp_problems(temp,solution.D,solution.Func_num);
	FVorigin=solution.Y[R[0]];
	if (solution.CECorBBOB==0)
		FVresult=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
	else if (solution.CECorBBOB==1)
		FVresult=bbob09.FunctionCalculation(X,solution.Func_num);
	else if (solution.CECorBBOB==2)
		FVresult=CEC13.cec13_problems(X,solution.D,solution.Func_num);
	else if (solution.CECorBBOB==3)
		FVresult=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
	
	//if using OBL, QOBL and QI, original initialization solutions are already evaluated and
	//should not be count as 1 FE if used in selection

	if(FVorigin<FVresult)
	{
		//do nothing, keep S unchanged
	}
	else
	{
		for (int i=0;i<solution.D;i++)
		{
			solution.S[R[0]][i]=X[i];
		}
		solution.Y[R[0]]=FVresult;
	}

	delete []temp;
}

void DE::DEprocess(int loop)
{
	extern Solution solution;
	SolutionOperator SOperator;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	FEoutput feoutput;
	
	int FE=0;//to count the FE times
	//int *R=new int[1];//random selected solution candidate
	//R[0]=0;
	int R=0;
	double *X=new double[solution.D];
	//feoutput.Initial();
	int count=0;
	
	SOperator.UpdateY();
	//if (solution.Tech_num==1||solution.Tech_num==2||solution.Tech_num==6)
	//{
	//	FE=FE+solution.P;
	//}
	
	double r1,r2;
	
	while(FE<=solution.MaxFE)
	{
		if (solution.AdaptiveFlag==1)
		{
			r1=rand();
			if (r1<0.1)
			{
				DEp.F=0.1+0.9*rand();
			}
			r2=rand();
			if (r2<0.1)
			{
				DEp.CR=rand();
			}
		}
		
		DEmutation(X,&R);
		DEcrossover(X,&R);
		DEselection(X,&R);
		feoutput.UpdateBest(FE);

		//record aggregated convergence
		if (FE%solution.NumConvergenceStep == 0)
		{
			solution.AvgConvergence[(FE/solution.NumConvergenceStep)] += solution.BestEver;
		}
		FE=FE+1;
	}
	
	feoutput.Release();
	delete []X;
	//delete []R;
}