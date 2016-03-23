//main function, designed for minimization problem
#include "stdafx.h"
/*#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include <windows.h>
#include <fstream>
#include "initialization.h"
#include "solution.h"
#include "Algorithm.h"
#include "FEoutput.h"
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "bbob09functions.h"*/
#include "bbob09supportfunctions.h"
//#include <vld.h>

using std::string;
using namespace std;

int main()
{
	srand(time(NULL));
	
	extern Solution solution;
	extern FEvariable fevariable;
	Initialization initialization;
	Algorithm algorithm;
	FEoutput feoutput;
	SolutionOperator SOperator;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	BBOB09support bbob09support;

	solution.Tech_num=1;
	solution.dir_file="new-joe-kuo-6.21201.txt";
	//*1 for random, 2 for chaos, 3 for opposition, 
	//*4 for Quasi Opposition, 5 for Quasi Interpolation, 6 for sobol
	//int Max_Tech=6;
	
	//int N=solution.Func_num*10; //dimension of each candidate
	//solution.MaxFE=500*solution.Func_num;//if =0, the result is the initialization result
	solution.Algorithm=1;
	//*1 for DE, 2 for CRO, 3 for PSO
	
	solution.ConstraintHandle="Bounceback";//Equaltobound/Bounceback
	solution.AdaptiveFlag=0;//1 turn on adaptive scheme, 0 turn off adaptive scheme

	//note if DE is used, P=5*N is suggested, revised inside function loop
	solution.TargetBest=0;
	solution.Error=0.00000001;
	solution.Looptime=1000;
	
	solution.trialID=1;  
	int StartFunction=1;//the first function that will be tested
	int MaxFunctionTested=1;//how many function is not being tested, 1 stands for 1 function.
	//Control variables: solution.D, MaxFE/D, RA
	solution.D=2;
	int IfUseBestP = 0;//if to use data in 2-D Array BestP, 0 is not to use, 1 is to use.
	int FixedEANPsize = 10; //less than 10 will trigger problem
	solution.EANP=FixedEANPsize;//reset P, CRO will change the population size
	//solution.seed = 1.0;//initiate random seed
	//solution.NumSeed=1;
	solution.NumMove=solution.EANP;
	solution.NeighborMethod = "Square";//to select the neighbor generation method, gaussian or square
	solution.NeighborCount = 10;//number of neighborhood tested
	solution.MinNeighborRange = 1000;//(Dmax-Dmin)/1000 -> minimum neighbor range
	solution.MaxNeighborRange = 10; //(Dmax-Dmin)/10 -> maximum neighbor range

	solution.CECorBBOB=2;//select which benchmark to use
	//0 for CEC14 expensive, 1 for BBOB09, 2 for CEC13, 3 for CEC14
	int RAorP=0;//select whether to use P or RA as x-axis
	//0 for P
	//1 for RA
	int FreeIni = 0;//decide if iniT is free or not
	//0 means non-free IniT
	//1 means free IniT, traditional assumption
	solution.CRODecomS = 5;

	//int FEarray[8] = {25,50,75,100,150,200,300,500};
	int FEarray[1] = {50};
	int MaxFEperDcount = 1;
	solution.IniNP=solution.EANP;//number of initial solutions
	solution.MaxFE=FEarray[0]*solution.D-solution.IniNP;
	
	solution.NumConvergenceStep = 10;
	solution.NumConvergence = FEarray[0]*solution.D/solution.NumConvergenceStep;
	bbob09.BBOBparameterInitial();

	//simulation information:
	cout<<"Technique: "<<solution.Tech_num<<endl;
	cout<<"Algorithm: "<<solution.Algorithm<<endl;
	cout<<"Function set: "<<solution.CECorBBOB<<endl;
	cout<<"Function tested: "<<StartFunction<<"-"<<StartFunction+MaxFunctionTested-1<<endl;
	cout<<"Loop time: "<<solution.Looptime<<endl;
	cout<<"BestP? "<<IfUseBestP;
	if (IfUseBestP==0)
		cout<<" Fixed NP: "<<FixedEANPsize<<endl;
	else
		cout<<endl;
	cout<<"RAorP? "<<RAorP<<endl;
	cout<<"FE tested: ";
	for (int i=0;i<MaxFEperDcount;i++)
	{
		cout<<FEarray[i]<<" ";
	}
	cout<<endl;
	
	solution.Min=new double[solution.D];
	solution.Max=new double[solution.D];
	solution.GlobalOptima=new double[solution.D];
	for (int d=0;d<solution.D;d++)
	{
		solution.GlobalOptima[d] = 0.0;
	}
	solution.LoopResults=new double[solution.Looptime];
	
	/*double ***IniRecord=new double**[solution.NumMove+1];//to record the averaged best initial solutions
	for (int i=0;i<solution.NumMove+1;i++)
	{
		IniRecord[i]=new double*[solution.EANP];
		for (int j=0;j<solution.EANP;j++)
		{
			IniRecord[i][j] = new double[solution.D];
		}
	}*/

	/*LARGE_INTEGER Frequency;//计数器频率  
	LARGE_INTEGER start_PerformanceCount;//起始计数器  
	LARGE_INTEGER end_PerformanceCount;//结束计数器 
	double run_time=0; //运行时间  
	QueryPerformanceFrequency(&Frequency);
	QueryPerformanceCounter(&start_PerformanceCount);*/
	int FunctionCount=0;

	int TotalProgress = solution.Looptime * MaxFunctionTested * solution.NeighborCount * (solution.NumMove+1);

	while(FunctionCount<MaxFunctionTested)//test one function first
	{
		//update function variables according to current Function count
		
		SOperator.UpdateFunction(StartFunction, FunctionCount, FixedEANPsize);
		
		FILE* loopfout_main;
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
		/*sprintf_s(FinalFileName,"ResultOutput//%s_A%d_Tech%d_Function%d_D%d_NumMove%d_%s.csv",FSet,solution.Algorithm,solution.Tech_num,solution.Func_num,solution.D,solution.NumMove,NMethod);
		fstream finalfoutclear(FinalFileName,ios::out);
		finalfoutclear.close();
		errno_t err=fopen_s(&loopfout_main,FinalFileName,"a+");
		if (err==0)
		{
			fprintf(loopfout_main,"Algorithm, %d,Tech, %d,D, %d\n Num of Moved, %d, Loop time, %d \n\n",solution.Algorithm,solution.Tech_num,solution.D, solution.NumMove, solution.Looptime);
		}
		fclose(loopfout_main);*/
		//define the range of "neighborhood". Currently only the min value is used 

		int neighbor_count = 0;

		//----------------------------------------------
		//start the seed loop
		while (neighbor_count<solution.NeighborCount+1)
		{
			//set the random seed, is this proper?
			//srand(log(seed_count*seed_count));
			srand(time(NULL));
			
			double NeighborhoodStep = (log10(solution.MinNeighborRange) - log10(solution.MaxNeighborRange))/solution.NeighborCount;
			double CurrentNeighborRange = log10(solution.MaxNeighborRange) + NeighborhoodStep*neighbor_count;
			
			for (int d=0;d<solution.D;d++)
			{
				solution.RangeofNeighborhood[d] = (solution.Max[d]-solution.Min[d])/(pow(10,CurrentNeighborRange));
				//solution.MaxRangeNeighbor[d] = (solution.Max[d]-solution.Min[d])/50;
			}

			//initialization
			initialization.Initial();
			//need to test the randomness of initial solutions
			for (int i=0;i<solution.EANP;i++)
			{
				for (int d=0;d<solution.D;d++)
				{
					solution.UsedIniS[i][d] = solution.S[i][d];
					solution.OriginalIniS[i][d] = solution.S[i][d];//original record the original initial solutions, will be used to update each move
				}
			}

			double best_ini_sum=0;
			//----------------------------------------------
			//start the move loop
			int move_count = 0;
			while (move_count<solution.NumMove+1)
			{
				initialization.MoveIniSolution(move_count);//here UsedIniS is changed
				//to record the average convergence
				for (int i=0;i<solution.NumConvergence;i++)
				{
					solution.AvgConvergence[i]=0;
				}

				int loop = 0;
				while(loop<solution.Looptime)
				{
					//If in future loop test is constructed, remember to reset the parameters as a new loop
					//initial solutions regarding to the given initialization technique and problem function		

					for (int i=0;i<solution.EANP;i++)
					{
						for (int d=0;d<solution.D;d++)
						{
							//copy back the same initial solutions (maybe some is moved)
							solution.S[i][d] = solution.UsedIniS[i][d];
							//cout<<solution.S[i][d]<<" ";
						}
						//cout<<endl;
					}
					//getchar();
					solution.BestEver=10000000000;
				
					best_ini_sum+=solution.Y[0];
					algorithm.Optimization(loop);
					solution.LoopResults[loop]=solution.BestEver;
			
					loop++;
					
					//progress display:
					int CurrentProgress = FunctionCount*solution.NeighborCount*(solution.NumMove+1)*solution.Looptime
						+ neighbor_count*(solution.NumMove+1)*solution.Looptime + move_count*solution.Looptime + loop;
					double Percentage = double(CurrentProgress)/double(TotalProgress);
					printf("%.2lf% % \r",Percentage*100);
					
					//to reset solution.EANP since CRO might change its value
					solution.EANP = FixedEANPsize;
				};//end of optimization loop
				
				//to record the average convergence
				for (int i=0;i<solution.NumConvergence;i++)
				{
					//average out
					solution.AvgConvergence[i]=solution.AvgConvergence[i]/double(solution.Looptime);
				}
				
				//get loop best, worst, mean and std, output to .csv together with convergence and initial solution
				feoutput.Loopsort();
				feoutput.Loopmeanstd(move_count, neighbor_count);
				
				move_count++;
			}
			neighbor_count++;
		}
		//-------------------------------------------

		feoutput.LOutput();
		feoutput.LOutputForR();

		if (solution.CECorBBOB==0)
		{
			solution.Func_num=solution.Func_num+3;
		}
		else
		{
			solution.Func_num++;
		}
		FunctionCount++;
		//cout<<endl;
		SOperator.Release(FixedEANPsize);
	}//end of problem loop

	//QueryPerformanceCounter(&end_PerformanceCount);
	//run_time = ( end_PerformanceCount.QuadPart - start_PerformanceCount.QuadPart ) / (double)Frequency.QuadPart;

	delete []solution.Min;
	delete []solution.Max;
	delete []solution.GlobalOptima;
	delete []solution.LoopResults;	
	bbob09.Release();
	/*for (int i=0;i<solution.NumMove;i++)
	{
		delete []IniRecord[i];
	}
	delete []IniRecord;*/

	cout<<"program end";
	getchar();
	return 0;
}