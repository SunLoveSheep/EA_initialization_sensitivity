#include "stdafx.h"
/*#include "solution.h"
#include <iostream>
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "bbob09functions.h"*/

using namespace std;

Solution solution;

SolutionOperator::SolutionOperator()
{
	
}

SolutionOperator::~SolutionOperator()
{
	
}

void SolutionOperator::UpdateFunction(int start_func, int func_count, int fixedEANPsize)
{
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	solution.EANP=fixedEANPsize;//EA population size
	solution.Pdouble=solution.CRODecomS*fixedEANPsize;//upper limit of population size. only initial 2P positions
	solution.S=new double*[solution.CRODecomS*fixedEANPsize];
	for (int p=0;p<solution.CRODecomS*fixedEANPsize;p++)
	{
		solution.S[p]=new double[solution.D];
		for (int d=0;d<solution.D;d++)
		{
			solution.S[p][d]=0;
		}
	}
	solution.IniS=new double*[solution.IniNP];
	for (int p=0;p<solution.IniNP;p++)
	{
		solution.IniS[p]=new double[solution.D];
		for (int d=0;d<solution.D;d++)
		{
			solution.IniS[p][d]=0;
		}
	}
	solution.UsedIniS=new double*[solution.EANP];
	solution.OriginalIniS=new double*[solution.EANP];
	for (int i=0;i<solution.EANP;i++)
	{
		solution.UsedIniS[i]=new double[solution.D];
		solution.OriginalIniS[i]=new double[solution.D];
		for (int d=0;d<solution.D;d++)
		{
			solution.UsedIniS[i][d]=0;
			solution.OriginalIniS[i][d]=0;
		}
	}
	solution.Y=new double[solution.CRODecomS*fixedEANPsize];
	solution.IniY=new double[solution.IniNP];
	solution.AvgConvergence=new double[solution.NumConvergence];
	//solution.MinRangeNeighbor = new double[solution.D];
	//solution.MaxRangeNeighbor = new double[solution.D];
	solution.RangeofNeighborhood = new double[solution.D];

	solution.MeanperMoveNeighbor = new double*[solution.NumMove+1];
	solution.StdperMoveNeighbor = new double*[solution.NumMove+1];
	for (int m=0;m<solution.NumMove+1;m++)
	{
		solution.MeanperMoveNeighbor[m] = new double[solution.NeighborCount+1];
		solution.StdperMoveNeighbor[m] = new double[solution.NeighborCount+1];
		for (int n=0;n<solution.NeighborCount+1;n++)
		{
			solution.MeanperMoveNeighbor[m][n] = 0;
			solution.StdperMoveNeighbor[m][n] = 0;
		}
	}

	//int loop=0;
	for (int i=0;i<solution.Looptime;i++)
	{
		solution.LoopResults[i]=0;
	}

	switch(solution.CECorBBOB)
	{
	case 0:
		if (solution.D==10)
		{
			solution.Func_num=start_func*3-2+func_count;//1 as sphere, etc.
		}
		else if (solution.D==30)
		{
			solution.Func_num=start_func*3-1+func_count;
		}
		else if (solution.D==100)
		{
			solution.Func_num=start_func*3+func_count;
		}

		//FileIO, read the M rotation matrix and o shift vector information
		CEC14.FileIO(solution.D,solution.Func_num);
		//different min max constraints for Function 13-18
		if((solution.Func_num>12)&&(solution.Func_num<16))
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.Min[n]=-32;
				solution.Max[n]=32;
			}
		}
		else if((solution.Func_num>15)&&(solution.Func_num<19))
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.Min[n]=-600;
				solution.Max[n]=600;
			}
		}
		else
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.Min[n]=-20;
				solution.Max[n]=20;
			}
		}
	case 1://initial setting for BBOB09
		solution.Func_num=start_func+func_count;
		bbob09.CalculateParameter(solution.Func_num);
		for (int n=0;n<solution.D;n++)
		{
			solution.Min[n]=-5;
			solution.Max[n]=5;
		}
	case 2://initial setting for CEC13
		solution.Func_num=start_func+func_count;
		CEC13.FileIO(solution.D,solution.Func_num);
		for (int n=0;n<solution.D;n++)
		{
			solution.Min[n]=-100;
			solution.Max[n]=100;
		}
	case 3://initial setting for CEC13
		solution.Func_num=start_func+func_count;
		CEC14_normal.FileIO(solution.D,solution.Func_num);
		for (int n=0;n<solution.D;n++)
		{
			solution.Min[n]=-100;
			solution.Max[n]=100;
		}
	}
}

void SolutionOperator::Release(int fixedEAPsize)
{
	delete []solution.Y;
	delete []solution.IniY;
	for (int p=0;p<solution.CRODecomS*fixedEAPsize;p++)
	{
		delete []solution.S[p];
	}
	for (int p=0;p<solution.IniNP;p++)
	{
		delete []solution.IniS[p];
	}
	delete []solution.IniS;
	delete []solution.S;
	delete []solution.UsedIniS;
	delete []solution.AvgConvergence;
	//delete []solution.MaxRangeNeighbor;
	//delete []solution.MinRangeNeighbor;
	delete []solution.RangeofNeighborhood;
	for (int m=0;m<solution.NumMove;m++)
	{
		delete []solution.MeanperMoveNeighbor[m];
		delete []solution.StdperMoveNeighbor[m];
	}
	delete []solution.MeanperMoveNeighbor;
	delete []solution.StdperMoveNeighbor;
}

void SolutionOperator::UpdateY()
{
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	//solution.Y=new double[solution.P];
	double *temp=new double[solution.D];

	for (int p=0;p<solution.EANP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			temp[n]=solution.S[p][n];
		}
				
		if (solution.CECorBBOB==0)//use CEC14 expensive
		{
			solution.Y[p]=CEC14.cec14_eotp_problems(temp,solution.D,solution.Func_num);
		}
		else if (solution.CECorBBOB==1) //use BBOB
		{
			solution.Y[p]=bbob09.FunctionCalculation(temp,solution.Func_num);
		}
		else if (solution.CECorBBOB==2) //use CEC13
		{
			solution.Y[p]=CEC13.cec13_problems(temp,solution.D,solution.Func_num);
		}
		else if (solution.CECorBBOB==3) //use CEC14
		{
			solution.Y[p]=CEC14_normal.cec14_problems(temp,solution.D,solution.Func_num);
		}
	}

	delete []temp;
	//delete temp;
}