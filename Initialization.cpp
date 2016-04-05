#include "stdafx.h"
#include <algorithm>
/*#include <string>
#include "solution.h"
#include "Initialization.h"
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "bbob09functions.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>*/

using std::string;
using namespace std;

Initialization::Initialization()
{

}

Initialization::~Initialization()
{
	
}

double Initialization::RandomGen(double min, double max)
{
	double Min = (min*1000000);
    double Max = (max*1000000);
    double Rand = rand()*rand();
    double Result;

	if (min!=max)
	{
		Result = ((int)(Rand)%(int)(Max-Min))+Min;
	}
	else
	{
		Result=min;
	}
    return Result/1000000.0;
}

double Initialization::GaussRandomGen(double miu, double sigma)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
         
    phase = 1 - phase;
	
	X=X*sigma+miu;
    return X;
}

//this function move initial solution into the neighborhood of the optima
void Initialization::MoveIniSolution(int movecount)
{
	extern Solution solution;
	//movecount can be taken as the number of initial solution that would be moved into the neighborhood
	//first, generate movecount nearby optima solutions
	double **tempIni = new double*[movecount+1];
	
	if (solution.NeighborMethod == "Gaussian")
	{
		for (int i=0;i<movecount+1;i++)
		{
			tempIni[i]=new double[solution.D];
			for (int d=0;d<solution.D;d++)
			{
				//a gaussian noise around the optima
				tempIni[i][d] = GaussRandomGen(solution.GlobalOptima[d],solution.RangeofNeighborhood[d]);
			}
		}
	}
	else if (solution.NeighborMethod == "Square")
	{
		for (int i=0;i<movecount+1;i++)
		{
			tempIni[i]=new double[solution.D];
			for (int d=0;d<solution.D;d++)
			{
				//a hyper-cube noise centered at the optima
				tempIni[i][d] = RandomGen(solution.GlobalOptima[d]-solution.RangeofNeighborhood[d],solution.GlobalOptima[d]+solution.RangeofNeighborhood[d]);
			}
		}
	}

	//then, replace movecount original initial solutions
	//1. Random replacement
	int *tempSeq = new int[solution.EANP];//to swap the sequence and random select movecount+1 solutions to be replaced
	int swap_count = 0;
	while (swap_count<solution.EANP)
	{
		int n1= (int)RandomGen(0,solution.EANP);
		int n2= (int)RandomGen(0,solution.EANP);
		
		for (int d=0;d<solution.D;d++)
		{
			double temp = solution.OriginalIniS[n1][d];
			solution.OriginalIniS[n1][d] = solution.OriginalIniS[n2][d];
			solution.OriginalIniS[n2][d] = temp;
		}
		//getchar();
		swap_count++;
	}

	//copy the swapped original solutions to the solution.UsedIniS:
	for (int i=0;i<solution.EANP;i++)
	{
		for (int d=0;d<solution.D;d++)
		{
			solution.UsedIniS[i][d] = solution.OriginalIniS[i][d];
		}
	}

	//replace the first movecount+1 initial solutions:
	for (int i=0;i<movecount;i++)
	{
		for (int d=0;d<solution.D;d++)
		{
			solution.UsedIniS[i][d] = tempIni[i][d];
		}
	}

	for (int i=0;i<movecount+1;i++)
	{
		delete []tempIni[i];
	}
	delete []tempIni;
	delete []tempSeq;
}

void MatchfromIniSolutionsToEASolutions(double *temparray)
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	double *tempX = new double[solution.D];
	double min=1000000000000000000;
	
	int *tempCount = new int[solution.IniNP];

	for (int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			tempX[n]=solution.IniS[p][n];
		}
		solution.IniY[p] = temparray[p];
	}
	
	sort(temparray,temparray+solution.IniNP);

	for (int p=0;p<solution.EANP;p++)
	{
		solution.Y[p]=temparray[p];
	}
	for (int p=0;p<solution.EANP;p++)
	{
		for (int q=0;q<solution.IniNP;q++)
		{
			if (temparray[p]==solution.IniY[q])
				tempCount[p]=q;
		}

		for (int d=0;d<solution.D;d++)
		{
			solution.S[p][d]=solution.IniS[tempCount[p]][d];
		}
	}

	delete []tempX;
	delete []tempCount;
}

void Initialization::CornerPointBasedSearch()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	//1. read corner points from the curves
	//initial the corner points. For each point, the two attributes are the range of neighborhood and the required number of points to bring improvement.
	double CornerPoints_D10_1Percent[5][2] = {{1.327,7},{1.435,4},{1.543,3},{1.812,2},{1.974,1}};
	int BestPoint = 2; // the point that best utilize each point to cover more spaces

	//2. Randomly generate EANP initial solutions
	double **tempIni = new double*[solution.EANP];
	double **tempIninotsorted = new double*[solution.EANP];
	double *tempY = new double[solution.EANP];
	double *tempX = new double[solution.D];
	double *tempYnotsorted = new double[solution.EANP];
	for(int p=0;p<solution.EANP;p++) //repeat for all solutions
	{
		tempIni[p] = new double[solution.D];
		tempIninotsorted[p] = new double[solution.D];
		for (int d=0;d<solution.D;d++)
		{
			tempIni[p][d] = RandomGen(solution.Min[d],solution.Max[d]);
			tempIninotsorted[p][d] = tempIni[p][d];
			tempX[d] = tempIni[p][d];
		}

		switch(solution.CECorBBOB)
		{
		case 0:
			tempY[p]=CEC14.cec14_eotp_problems(tempX,solution.D,solution.Func_num);
			break;
		case 1:
			tempY[p]=bbob09.FunctionCalculation(tempX,solution.Func_num);
			break;
		case 2:
			tempY[p]=CEC13.cec13_problems(tempX,solution.D,solution.Func_num);
			break;
		case 3:
			tempY[p]=CEC14_normal.cec14_problems(tempX,solution.D,solution.Func_num);
			break;
		}
		tempYnotsorted[p] = tempY[p];
	}

	//3. sort the solutions
	sort(tempY,tempY+solution.EANP);
	//the below sorting algorithm cannot handle when there are more than two objective values are identical (though this is very unlikely the situation)
	for (int p=0;p<solution.EANP;p++)
	{
		for (int i=0;i<solution.EANP;i++)
		{
			if (tempYnotsorted[p] == tempY[i])
			{
				for (int d=0;d<solution.D;d++)
				{
					tempIni[i][d] = tempIninotsorted[p][d];
				}
				break;
			}
		}
	}
	//4. while IniFE is still available:
	int FECount = 0;
	double **NewPoints = new double*[(int)CornerPoints_D10_1Percent[BestPoint][1]-1];
	double **NewPointsExploration = new double*[solution.EANP-1];
	for (int i = 0;i<(int)CornerPoints_D10_1Percent[BestPoint][1]-1;i++)
	{
		NewPoints[i] = new double[solution.D];
	}
	for (int i = 0;i<solution.EANP-1;i++)
	{
		NewPointsExploration[i] = new double[solution.D];
	}
	while (FECount<solution.IniNP)
	{
		//4.1 do neighbor search on the best initial solution
		double NeighborRange = 0;
		double tempPointY = 0;
		for (int i = 0;i<(int)CornerPoints_D10_1Percent[BestPoint][1]-1;i++)
		{
			for (int d=0;d<solution.D;d++)
			{
				NeighborRange = (solution.Max[d] - solution.Min[d])/pow(10,CornerPoints_D10_1Percent[BestPoint][0]);
				NewPoints[i][d] = RandomGen(tempIni[0][d]-NeighborRange,tempIni[0][d]+NeighborRange);
				tempX[d] = NewPoints[i][d];
			}

			switch(solution.CECorBBOB)
			{
			case 0:
				tempPointY=CEC14.cec14_eotp_problems(tempX,solution.D,solution.Func_num);
				break;
			case 1:
				tempPointY=bbob09.FunctionCalculation(tempX,solution.Func_num);
				break;
			case 2:
				tempPointY=CEC13.cec13_problems(tempX,solution.D,solution.Func_num);
				break;
			case 3:
				tempPointY=CEC14_normal.cec14_problems(tempX,solution.D,solution.Func_num);
				break;
			}
			FECount++;//update FE count

			//if better new point is found, update
			if (tempPointY < tempY[0])
			{
				tempY[0] = tempPointY;
				for (int d=0;d<solution.D;d++)
				{
					tempIni[0][d] = NewPoints[i][d];
				}
			}
		}

		//4.2 do exploration on the other initial solutions
		for (int i=0;i<solution.EANP-1;i++)//start from 1, the first (smallest) has been used for neighbor search.
		{
			for (int d=0;d<solution.D;d++)//do exploration on each dimension. If lower(upper) bound of neighbor range is smaller than min[d](max[d]), then only random generat in the upper(lower) side
			//if both sides feasible, random generate will consider both sides
			//if either sides infeasible, return error message: neighborhood too big, no feasible solution.
			{
				double Lower, Upper;
				NeighborRange = (solution.Max[d] - solution.Min[d])/pow(10,CornerPoints_D10_1Percent[BestPoint][0]);
				Lower = tempIni[i+1][d] - NeighborRange - solution.Min[d];
				Upper = solution.Max[d] - (tempIni[i+1][d] + NeighborRange);
				if (Lower<0 && Upper <0)//both sides infeasible
				{
					printf_s("Too big neighbor range, no feasible solution on Dimension %d",d);
				}
				else if (Lower>0 && Upper>0)//both sides feasible
				{
					double r = rand();
					if (r<0.5)//generate at lower part
						tempX[d] = RandomGen(solution.Min[d], tempIni[i+1][d] - NeighborRange);
					else
						tempX[d] = RandomGen(tempIni[i+1][d] + NeighborRange, solution.Max[d]);
				}
				else if(Lower>0 && Upper<0)//only at lower part feasible
					tempX[d] = RandomGen(solution.Min[d], tempIni[i+1][d] - NeighborRange);
				else if(Lower<0 && Upper>0)//only upper part feasible
					tempX[d] = RandomGen(tempIni[i+1][d] + NeighborRange, solution.Max[d]);
			}

			//calculate objective function value:
			switch(solution.CECorBBOB)
			{
			case 0:
				tempPointY=CEC14.cec14_eotp_problems(tempX,solution.D,solution.Func_num);
				break;
			case 1:
				tempPointY=bbob09.FunctionCalculation(tempX,solution.Func_num);
				break;
			case 2:
				tempPointY=CEC13.cec13_problems(tempX,solution.D,solution.Func_num);
				break;
			case 3:
				tempPointY=CEC14_normal.cec14_problems(tempX,solution.D,solution.Func_num);
				break;
			}
			FECount++;//update FE count

			//update solution:
			if (tempPointY<tempY[i+1])
			{
				for (int d=0;d<solution.D;d++)
				{
					tempIni[i+1][d] = tempX[d];
				}
				tempY[i+1] = tempPointY;
			}
		}
		//question 1: how to balance the rate of these exploitation and exploration?
		//question 2: how to select new points? Greedy?

		//4.3 after exploitation and exploration, sort the solutions again:
		sort(tempY,tempY+solution.EANP);
		//the below sorting algorithm cannot handle when there are more than two objective values are identical (though this is very unlikely the situation)
		for (int p=0;p<solution.EANP;p++)
		{
			for (int i=0;i<solution.EANP;i++)
			{
				if (tempYnotsorted[p] == tempY[i])
				{
					for (int d=0;d<solution.D;d++)
					{
						tempIni[i][d] = tempIninotsorted[p][d];
					}
					break;
				}
			}
		}
	}
	
	//5. Copy the final initial solutions to the solution.Y and solution.S arrays.

	delete []tempY;
	delete []tempYnotsorted;
	delete []tempX;
	for(int p=0;p<solution.EANP;p++) //repeat for all solutions
	{
		delete []tempIni[p];
		delete []tempIninotsorted[p];
	}
	delete []tempIni;
	delete []tempIninotsorted;
	for (int i = 0;i<(int)CornerPoints_D10_1Percent[BestPoint][1]-1;i++)
	{
		delete []NewPoints[i];
	}
	for (int i = 0;i<solution.EANP-1;i++)
	{
		delete []NewPointsExploration[i];
	}
	delete []NewPoints;
	delete []NewPointsExploration;
}

//this function returns a randomly swapped integer vector containing values from Min to Max. Min should be 0 and Max is data.NumSample-1
void RandomIntVec(int Min, int Max, int *Vec)
{
	extern Solution solution;
	for (int i=0;i<solution.IniNP;i++){
		Vec[i]=i;
	}

	//swap two random numbers in the vec for NumSample times
	int swap_time = 0;
	int swap_temp = 0;
	while (swap_time<solution.IniNP){
		int rNum1 = (int)(rand()*rand())%(Max+1-Min)+Min-1;//random generate two integers
		int rNum2 = (int)(rand()*rand())%(Max+1-Min)+Min-1;//random generate two integers
		//swap the numbers on that integer position of Vec
		swap_temp = Vec[rNum1];
		Vec[rNum1] = Vec[rNum2];
		Vec[rNum2] = swap_temp;

		swap_time++;//counter update
	}
}

//this function generate index samples for LHD
void LHD_IndexSamples(int **OutIndex)
{
	extern Solution solution;
	
	int *temp_vec = new int[solution.IniNP];
	
	for (int d=0;d<solution.D;d++)
	{
		RandomIntVec(1,solution.IniNP,temp_vec);
		for (int n=0;n<solution.IniNP;n++)
		{
			OutIndex[n][d] = temp_vec[n];
		}
	}

	delete []temp_vec;
}

void LHD_RealSampling(int **InIndex, double **OutSamples)//this function takes index samples and generate actual samples by means of LHD
{
	extern Solution solution;
	Initialization initialization;
	
	for (int d=0;d<solution.D;d++)
	{
		double Sec_Length = (double)((solution.Max[d]-solution.Min[d]))/(double)(solution.IniNP);//to calculate the length of each section, dividing according to the number of samples
		
		for (int n=0;n<solution.IniNP;n++)
		{
			//randomly generate real sample values according to index.
			double randnum = initialization.RandomGen(0.0,Sec_Length);
			OutSamples[n][d] = solution.Min[d] + (InIndex[n][d]) * Sec_Length +randnum;
		}
	}
}

void Initialization::LHD()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	//initiate variables
	int **LHD_index = new int*[solution.IniNP];
	double **LHD_samples = new double*[solution.IniNP];
	for (int i=0;i<solution.IniNP;i++)
	{
		LHD_index[i] = new int[solution.D];
		LHD_samples[i] = new double[solution.D];
		for (int d=0;d<solution.D;d++)
		{
			LHD_index[i][d]=0;
			LHD_samples[i][d]=0;
		}
	}

	//LHD indexing and sampling
	LHD_IndexSamples(LHD_index);
	LHD_RealSampling(LHD_index,LHD_samples);
	/*if (data.SequencingMethod == "NN")
	{
		//LHD samples are randomly generated, use nearest neighbor to sort and sequencing samples
		NearestSorting(LHD_samples, data.SampleSequence);
	}*/
	//else if (data.SequencingMethod == "RNG")
	//{
		//LHD samples are already randomly generated. Directly copy to output sequence
	for (int i=0;i<solution.IniNP;i++)
	{
		for (int d=0;d<solution.D;d++)
		{
			solution.IniS[i][d] = LHD_samples[i][d];
		}
	}
	//}
	
	double *X=new double[solution.D];
	double *Y=new double[solution.IniNP];
	for(int p=0;p<solution.IniNP;p++) //repeat for all solutions
	{
		for (int i=0;i<solution.D;i++)
		{
			X[i]=solution.IniS[p][i];
		}
		
		switch (solution.CECorBBOB)
		{
		case 0:
			Y[p]=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
			break;
		case 1:
			Y[p]=bbob09.FunctionCalculation(X,solution.Func_num);
			break;
		case 2:
			Y[p]=CEC13.cec13_problems(X,solution.D,solution.Func_num);
			break;
		case 3:
			Y[p]=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
			break;
		}
	}

	MatchfromIniSolutionsToEASolutions(Y);
	
	//release RAM
	for (int i=0;i<solution.IniNP;i++)
	{
		delete []LHD_index[i];
		delete []LHD_samples[i];
	}
	delete []LHD_index;
	delete []LHD_samples;
	delete []X;
	delete []Y;
}

void Initialization::Random()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	double *X=new double[solution.D];
	double *Y=new double[solution.IniNP];
	for(int p=0;p<solution.IniNP;p++) //repeat for all solutions
	{
		for (int i=0;i<solution.D;i++)
		{
			solution.IniS[p][i]=RandomGen(solution.Min[i],solution.Max[i]);
			X[i]=solution.IniS[p][i];
		}

		switch(solution.CECorBBOB)
		{
		case 0:
			Y[p]=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
			break;
		case 1:
			Y[p]=bbob09.FunctionCalculation(X,solution.Func_num);
			break;
		case 2:
			Y[p]=CEC13.cec13_problems(X,solution.D,solution.Func_num);
			break;
		case 3:
			Y[p]=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
			break;
		}
	}

	MatchfromIniSolutionsToEASolutions(Y);
	delete []X;
	delete []Y;
}

void Initialization::Chaos()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	double *X=new double[solution.D];
	double *Y=new double[solution.IniNP];
	
	//numerical precision fixed at .xxxxxxx
	cout << setiosflags(ios::fixed) << setprecision(7);

	double **R=new double*[solution.IniNP];
	for (int p=0;p<solution.IniNP;p++)
	{
		R[p]=new double[solution.D];
	}

	for(int p=0;p<solution.IniNP;p++) //repeat for all solutions
	{
		R[p][0]=RandomGen(0,1); //for each solution, the first variable is randomly initialized
		for (int i=1;i<solution.D;i++)
		{
			R[p][i]=4*R[p][i-1]*(1-R[p][i-1]);//for the rest variables of the solution, chaosly assign values
		}
	}

	//to recover R from [0,1] to the real range [min,max]
	for(int p=0;p<solution.IniNP;p++) //repeat for all solutions
	{
		for (int i=0;i<solution.D;i++)
		{
			solution.IniS[p][i]=solution.Min[i]+(solution.Max[i]-solution.Min[i])*R[p][i];
			X[i]=solution.IniS[p][i];
		}
		if (solution.CECorBBOB==0)
			Y[p]=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			Y[p]=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			Y[p]=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			Y[p]=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
	}

	MatchfromIniSolutionsToEASolutions(Y);

	//return R;
	for (int p=0;p<solution.IniNP;p++)
	{
		delete []R[p];
	}
	delete []R;
	delete []X;
	delete []Y;
}

void Initialization::Opposition()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	double *X=new double[solution.D];
	double *Y=new double[solution.IniNP];

	double tempY,tempYop; //objective function value of original and its opposition solutions
	for(int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X[n]=RandomGen(solution.Min[0],solution.Max[0]);//randomly generate the variables
		}
		if (solution.CECorBBOB==0)
			tempY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempY=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempY=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempY=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);

		for (int n=0;n<solution.D;n++)
		{
			X[n]=solution.Max[n]-X[n]+solution.Min[n];//get the opposition of the solution
		}
		if (solution.CECorBBOB==0)
			tempYop=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempYop=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempYop=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempYop=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);

		//suppose minimization problem
		if(tempY<=tempYop)//if original better, keep the original solution
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.IniS[p][n]=solution.Max[n]-X[n]+solution.Min[n];
				//since solution.X already been opposited, opposite it back to get the original solution
			}
			Y[p]=tempY;
		}
		else// if opposition result better, use the opposition solution.
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.IniS[p][n]=X[n];
			}
			Y[p]=tempYop;
		}
	}

	MatchfromIniSolutionsToEASolutions(Y);

	delete []X;
	delete []Y;
}

void Initialization::QuasiOpposition()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	/*double **R=new double*[solution.Pdouble];//times 2, in order that the population size is changed in some algorithms
	for (int p=0;p<solution.Pdouble;p++)
	{
		R[p]=new double[solution.D];
	}*/
	
	double *X=new double[solution.D];
	double *QX=new double[solution.D];
	double *Y=new double[solution.IniNP];

	double Middle;//record the middle point

	double tempY,tempYop; //objective function value of original and its opposition solutions
	for(int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X[n]=RandomGen(solution.Min[0],solution.Max[0]);//randomly generate the variables
		}
		if (solution.CECorBBOB==0)
			tempY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempY=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempY=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempY=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);

		for (int n=0;n<solution.D;n++)
		{
			//get the opposition of the solution
			QX[n]=solution.Max[n]-X[n]+solution.Min[n];
			Middle=(solution.Min[n]+solution.Max[n])/2;

			//quasi scheme
			if(X[n]<Middle)
			{
				QX[n]=Middle+(QX[n]-Middle)*rand();
			}
			else
			{
				QX[n]=QX[n]+(Middle-QX[n])*rand();
			}	
		}
		if (solution.CECorBBOB==0)
			tempYop=CEC14.cec14_eotp_problems(QX,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempYop=bbob09.FunctionCalculation(QX,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempYop=CEC13.cec13_problems(QX,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempYop=CEC14_normal.cec14_problems(QX,solution.D,solution.Func_num);

		//suppose minimization problem
		if(tempY<=tempYop)//if original better, keep the original solution
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.S[p][n]=X[n];
				//since solution.X already been opposited, opposite it back to get the original solution
			}
			Y[p]=tempY;
		}
		else// if opposition result better, use the opposition solution.
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.S[p][n]=QX[n];
			}
			Y[p]=tempYop;
		}
	}

	MatchfromIniSolutionsToEASolutions(Y);

	delete []X;
	delete []QX;
	delete []Y;
}

void Initialization::QuasiInterpolation()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	int DoubleIniNP=2*solution.IniNP;
	double **R=new double*[DoubleIniNP];//times 2, in order that the population size is changed in some algorithms
	for (int p=0;p<DoubleIniNP;p++)
	{
		R[p]=new double[solution.D];
	}
	/*double **Rfinal=new double*[solution.Pdouble];//times 2, in order that the population size is changed in some algorithms
	for (int p=0;p<solution.Pdouble;p++)
	{
		Rfinal[p]=new double[solution.D];
	}*/
	
	double *X=new double[solution.D];
	double *Y=new double[DoubleIniNP];
	double *rY=new double[solution.IniNP];
	
	double tempY;
	double Min=10000000000000;
	int MinCount=0;

	for(int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X[n]=RandomGen(solution.Min[n],solution.Max[n]);//randomly generate the variables
			R[p][n]=X[n];
		}
		if (solution.CECorBBOB==0)
			tempY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			tempY=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			tempY=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			tempY=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
		solution.IniY[p]=tempY;
		
		//record the best random
		if(tempY<Min)
		{
			Min=tempY;
			MinCount=p;
		}
	}

	int b,c;
	double bY,cY,MinCountY;
	for (int p=solution.IniNP;p<DoubleIniNP;p++)
	{
		//find b and c different from MinCount
		b=rand()%(solution.IniNP+1);
		while(b==MinCount)
		{
			b=rand()%(solution.IniNP+1);
		}
		c=rand()%(solution.IniNP+1);
		while((c==b)||(c==MinCount))
		{
			c=rand()%(solution.IniNP+1);
		}
		//copy the objective function value
		for (int n=0;n<solution.D;n++)
		{
			X[n]=R[b][n];	
		}
		//bY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		bY=solution.IniY[b];
		for (int n=0;n<solution.D;n++)
		{
			X[n]=R[c][n];	
		}
		//cY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		cY=solution.IniY[c];
		for (int n=0;n<solution.D;n++)
		{
			X[n]=R[MinCount][n];	
		}
		//MinCountY=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		MinCountY=solution.IniY[MinCount];

		//calculate the quasi interpolation points
		for (int n=0;n<solution.D;n++)
		{
			R[p][n]=0.5*(((R[b][n]*R[b][n]-R[c][n]*R[c][n])*MinCountY+(R[c][n]*R[c][n]-R[MinCount][n]*R[MinCount][n])*bY+(R[MinCount][n]*R[MinCount][n]-R[b][n]*R[b][n])*cY)
				/((R[b][n]-R[c][n])*MinCountY+(R[c][n]-R[MinCount][n])*bY+(R[MinCount][n]-R[b][n])*cY));
		}	
	}
	
	//select the best solution.P solutions
	for (int p=0;p<solution.IniNP;p++)
	{
		Y[p]=solution.IniY[p];
	}
	for (int p=solution.IniNP;p<DoubleIniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X[n]=R[p][n];
		}
		//The following change can hugely reduce the RAM occupied
		if (solution.CECorBBOB==0)
			Y[p]=CEC14.cec14_eotp_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==1)
			Y[p]=bbob09.FunctionCalculation(X,solution.Func_num);
		else if (solution.CECorBBOB==2)
			Y[p]=CEC13.cec13_problems(X,solution.D,solution.Func_num);
		else if (solution.CECorBBOB==3)
			Y[p]=CEC14_normal.cec14_problems(X,solution.D,solution.Func_num);
	}

	int *Ycount=new int[DoubleIniNP];//to record p after sorting Y
	for (int i=0; i<solution.Pdouble; i++)
	{
		Ycount[i]=i;
	}
	//sorting Y
	for (int i=0; i<solution.Pdouble; i++)
	{  
        for (int j=solution.Pdouble-1; j>i; j--)  
        {  
			if (Y[j]<Y[j-1])
			{
				swap(Y[j], Y[j-1]);
				swap(Ycount[j],Ycount[j-1]);
			}
		}
	}

	for (int p=0;p<solution.IniNP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			solution.IniS[p][n]=R[Ycount[p]][n];
		}
		rY[p]=Y[p];
	}

	MatchfromIniSolutionsToEASolutions(rY);

	for (int p=0;p<solution.Pdouble;p++)
	{
		delete []R[p];
	}
	delete []R;
	/*for (int p=0;p<solution.Pdouble;p++)
	{
		delete []Rfinal[p];
	}
	delete []Rfinal;*/
	delete []X;
	delete []Y;
	delete []rY;
	delete []Ycount;
}

void Initialization::Sobol()
{
	extern Solution solution;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;

	int omit_N = 5; //omit the first omit_N numbers
	int generate_P = 0;
	generate_P = omit_N + 2*solution.IniNP;
	ifstream infile(solution.dir_file,ios::in);
	if (!infile) {
    cout << "Input file containing direction numbers cannot be found!\n";
	getchar();
    exit(1);
	}
	char buffer[1000];
  infile.getline(buffer,1000,'\n');
  
  // L = max number of bits needed 
  unsigned L = (unsigned)ceil(log((double)generate_P)/log(2.0)); 
  
  // C[i] = index from the right of the first zero bit of i
  unsigned *C = new unsigned [generate_P];
  C[0] = 1;
  for (int i=1;i<=generate_P-1;i++) {
    C[i] = 1;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // POINTS[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  double **POINTS = new double * [generate_P];
  for (int i=0;i<=generate_P-1;i++) POINTS[i] = new double [solution.D];
  for (int j=0;j<=solution.D-1;j++) POINTS[0][j] = 0; 
  
  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
  unsigned *V = new unsigned [L+1]; 
  for (unsigned i=1;i<=L;i++) V[i] = 1 << (32-i); // all m's = 1

  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  unsigned *X = new unsigned [generate_P];
  X[0] = 0;
  for (int i=1;i<=generate_P-1;i++) {
    X[i] = X[i-1] ^ V[C[i-1]];
    POINTS[i][0] = (double)X[i]/pow(2.0,32); // *** the actual points
    //        ^ 0 for first dimension
  }
  
  // Clean up
  delete [] V;
  delete [] X;
    
  // ----- Compute the remaining dimensions -----
  for (int j=1;j<=solution.D-1;j++) {
    
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i=1;i<=s;i++) infile >> m[i];
    
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    unsigned *V = new unsigned [L+1];
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i] = m[i] << (32-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i] = m[i] << (32-i); 
      for (unsigned i=s+1;i<=L;i++) {
	V[i] = V[i-s] ^ (V[i-s] >> s); 
	for (unsigned k=1;k<=s-1;k++) 
	  V[i] ^= (((a >> (s-1-k)) & 1) * V[i-k]); 
      }
    }
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [generate_P];
    X[0] = 0;
	
    for (int i=1;i<=generate_P-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
	  
      POINTS[i][j] = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
   }
    // Clean up
    delete [] m;
    delete [] V;
    delete [] X;
  }
  delete [] C;
    
  double *Y=new double[solution.IniNP];
  double *rX=new double[solution.D];
  for (int p=0;p<solution.Pdouble;p++)
  {
	for (int n=0;n<solution.D;n++)
	{
		solution.IniS[p][n]=POINTS[p+omit_N][n];
		rX[n]=solution.IniS[p][n];
	}
	if (solution.CECorBBOB==0)
		Y[p]=CEC14.cec14_eotp_problems(rX,solution.D,solution.Func_num);
	else if (solution.CECorBBOB==1)
		Y[p]=bbob09.FunctionCalculation(rX,solution.Func_num);
	else if (solution.CECorBBOB==2)
		Y[p]=CEC13.cec13_problems(rX,solution.D,solution.Func_num);
	else if (solution.CECorBBOB==3)
		Y[p]=CEC14_normal.cec14_problems(rX,solution.D,solution.Func_num);
  }

  MatchfromIniSolutionsToEASolutions(Y);

  for (int p=0;p<generate_P;p++)
  {
	delete []POINTS[p];
  }
  delete []POINTS;
}

void Initialization::Initial()
{
	extern Solution solution;

	switch(solution.Tech_num)
	{
	case 1:
		Random();
		break;
	case 2:
		Chaos();
		break;
	case 3:
		Opposition();
		break;
	case 4:
		QuasiOpposition();
		break;
	case 5:
		QuasiInterpolation();
		break;
	case 6:
		Sobol();
		break;
	case 7:
		LHD();
		break;
	}
}
