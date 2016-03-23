#include "stdafx.h"
//#include "solution.h"
//#include "FEoutput.h"
//#include <iostream>
//#include <math.h>
//#include "cec14_eotp.h"

using namespace std;

FEvariable fevariable;

FEoutput::FEoutput()
{
	extern Solution solution;
	
	Loopbest=0;
	Loopworst=0;
	Loopmean=0;
	Loopstd=0;
}

FEoutput::~FEoutput()
{
	
}

void FEoutput::Initial()
{
	extern Solution solution;

	//Bestever=new double[solution.D];
	BesteverFE=0;
}

void FEoutput::Release()
{
	//delete []Bestever;
}

void FEoutput::UpdateBest(int FE)
{
	extern Solution solution;
	
	for (int p=0;p<solution.EANP;p++)
	{
		if (solution.BestEver>solution.Y[p])
		{
			solution.BestEver=solution.Y[p];
			//for (int n=0;n<solution.D;n++)
			//{
			//	Bestever[n]=solution.S[p][n];
			//}
			//BesteverFE=FE;
		}
	}
}

void FEoutput::Loopsort()
{	
	extern Solution solution;
	
	double max=-10000000000000, min=1000000000000000;
	
	for(int i=0;i<solution.Looptime;i++)
	{
		if(solution.LoopResults[i]>max)
		{
			max=solution.LoopResults[i];
		}
		if(solution.LoopResults[i]<min)
		{
			min=solution.LoopResults[i];
		}
	}

	Loopbest=min;
	Loopworst=max;	
}

void FEoutput::Loopmeanstd(int movecount, int neighborcount)
{
	extern Solution solution;
	
	double sum=0;
	for (int i=0;i<solution.Looptime;i++)
	{
		sum=sum+solution.LoopResults[i];
	}
	Loopmean=sum/solution.Looptime;

	double sumstd=0;
	for (int i=0;i<solution.Looptime;i++)
	{
		sumstd=sumstd+(solution.LoopResults[i]-Loopmean)*(solution.LoopResults[i]-Loopmean);
	}
	Loopstd=sqrt(sumstd/(double)solution.Looptime);

	//cout<<endl<<sumstd<<" "<<solution.Looptime<<" "<<sumstd/(double)solution.Looptime<<" "<<sqrt(sumstd/(double)solution.Looptime)<<endl;
	//cout<<"loop mean: "<<Loopmean<<endl;
	//getchar();

	solution.MeanperMoveNeighbor[movecount][neighborcount] = Loopmean;
	solution.StdperMoveNeighbor[movecount][neighborcount] = Loopstd;

}

void FEoutput::LOutput()
{
	extern Solution solution;

	FILE* loopfout;
	char FinalFileName[100];
	const char *FSet;
	switch(solution.CECorBBOB)
	{
	case 0:
		FSet = "CEC14_eotp";
	case 1:
		FSet = "BBOB'09";
	case 2:
		FSet = "CEC13";
	case 3:
		FSet = "CEC14";
	}
	const char* NMethod = solution.NeighborMethod.c_str();
	sprintf_s(FinalFileName,"ResultOutput//%s_A%d_Tech%d_Function%d_D%d_NumMove%d_%s.csv",FSet,solution.Algorithm,solution.Tech_num,solution.Func_num,solution.D,solution.NumMove,NMethod);
	//fstream finalfoutclear(FinalFileName,ios::out);
	//finalfoutclear.close();
	errno_t err=fopen_s(&loopfout,FinalFileName,"a+");
	if (err!=0)
	{
		exit(1);
	}
	else
		fprintf(loopfout,"Algorithm, %d,Tech, %d,D, %d\n Num of Moved, %d, Loop time, %d \n\n",solution.Algorithm,solution.Tech_num,solution.D, solution.NumMove, solution.Looptime);


	for (int n=0;n<solution.NeighborCount+1;n++)
	{
		//calculate currrent neighborhood range and write to csv file
		double NeighborhoodStep = (log10(solution.MinNeighborRange) - log10(solution.MaxNeighborRange))/solution.NeighborCount;
		//double CurrentNeighborRange = pow(10,(log10(solution.MaxNeighborRange) + NeighborhoodStep*n));
		double CurrentNeighborRange = (log10(solution.MaxNeighborRange) + NeighborhoodStep*n);
		fprintf_s(loopfout,"NeighborRange:, %f \n Num of Moved:,",CurrentNeighborRange);

		//for each neighborhood range, output all number of moved results:
		for (int m=0;m<solution.NumMove+1;m++)
		{
			fprintf_s(loopfout,"%d,",m);
		}
		fprintf_s(loopfout,"\n Neighbor %f LoopMean,",CurrentNeighborRange);
		//print loop mean for each moved and each neighborhood range
		for (int m=0;m<solution.NumMove+1;m++)
		{
			fprintf_s(loopfout,"%f,",solution.MeanperMoveNeighbor[m][n]);
		}
		fprintf_s(loopfout,"\n Neighbor %f LoopStd,",CurrentNeighborRange);
		//print loop mean for each moved and each neighborhood range
		for (int m=0;m<solution.NumMove+1;m++)
		{
			fprintf_s(loopfout,"%f,",solution.StdperMoveNeighbor[m][n]);
		}
		fprintf_s(loopfout,"\n\n");
	}

	//print initial solution
	/*for (int i=0;i<solution.EANP;i++)
	{
		for (int d=0;d<solution.D;d++)
		{
			fprintf_s(loopfout, "%f,",solution.UsedIniS[i][d]);
		}
		fprintf_s(loopfout,"\n");
	}*/

	//print loop best, worst, mean, std
	//fprintf_s(loopfout,"\n Loop Best:, %f \n Loop Worst:, %f \n Loopmean:, %f \n Loopstd:, %f \n",Loopbest, Loopworst, Loopmean, Loopstd);
	//fprintf_s(loopfout,"\n Loopmean:, %f \n Loopstd:, %f \n", Loopmean, Loopstd);

	//print convergence:
	/*fprintf_s(loopfout, "Convergence:,\n FE,");
	for (int i=0;i<solution.NumConvergence;i++)
	{
		fprintf_s(loopfout,"%d,",i*solution.NumConvergenceStep);
	}
	fprintf_s(loopfout,"\n BestEver,");
	for (int i=0;i<solution.NumConvergence;i++)
	{
		fprintf_s(loopfout,"%f,",solution.AvgConvergence[i]);
	}
	fprintf_s(loopfout,"\n\n");*/

	if (loopfout == NULL)
	{
		cout<<"file is NULL!";
		getchar();
	}
	fclose(loopfout);
}

void FEoutput::LOutputForR()
{
	extern Solution solution;

	FILE* loopfoutforR;
	char FinalFileName[100];
	const char *FSet;
	switch(solution.CECorBBOB)
	{
	case 0:
		FSet = "CEC14_eotp";
	case 1:
		FSet = "BBOB'09";
	case 2:
		FSet = "CEC13";
	case 3:
		FSet = "CEC14";
	}
	const char* NMethod = solution.NeighborMethod.c_str();
	sprintf_s(FinalFileName,"ResultOutput//ForR_%s_A%d_Tech%d_Function%d_D%d_NumMove%d_%s.csv",FSet,solution.Algorithm,solution.Tech_num,solution.Func_num,solution.D,solution.NumMove,NMethod);
	fstream finalfoutclear(FinalFileName,ios::out);
	finalfoutclear.close();
	errno_t err=fopen_s(&loopfoutforR,FinalFileName,"a+");
	if (err!=0)
	{
		exit(1);
	}
	fprintf_s(loopfoutforR,"Num of Moved,");
		//for each neighborhood range, output all number of moved results:
	for (int m=0;m<solution.NumMove+1;m++)
	{
		fprintf_s(loopfoutforR,"%d,",m);
	}
	for (int n=0;n<solution.NeighborCount+1;n++)
	{
		//calculate currrent neighborhood range and write to csv file
		double NeighborhoodStep = (log10(solution.MinNeighborRange) - log10(solution.MaxNeighborRange))/solution.NeighborCount;
		//double CurrentNeighborRange = pow(10,(log10(solution.MaxNeighborRange) + NeighborhoodStep*n));
		double CurrentNeighborRange = (log10(solution.MaxNeighborRange) + NeighborhoodStep*n);
		//fprintf_s(loopfoutforR,"\n NeighborRange, %f",CurrentNeighborRange);
		fprintf_s(loopfoutforR,"\n Neighbor %f LoopMean,",CurrentNeighborRange);
		//print loop mean for each moved and each neighborhood range
		for (int m=0;m<solution.NumMove+1;m++)
		{
			fprintf_s(loopfoutforR,"%f,",solution.MeanperMoveNeighbor[m][n]);
		}
		fprintf_s(loopfoutforR,"\n Neighbor %f LoopStd,",CurrentNeighborRange);
		//print loop mean for each moved and each neighborhood range
		for (int m=0;m<solution.NumMove+1;m++)
		{
			fprintf_s(loopfoutforR,"%f,",solution.StdperMoveNeighbor[m][n]);
		}
		//fprintf_s(loopfoutforR,"\n");
	}

	if (loopfoutforR == NULL)
	{
		cout<<"file is NULL!";
		getchar();
	}
	fclose(loopfoutforR);
}

/*void FEoutput::LoopOprocess()
{
	Loopsort();
	Loopmeanstd();
	//LOutput();

	//delete []Loopresults;
}*/