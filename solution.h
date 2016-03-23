//solution candidate definition
#include <string>
#include <vector>

using std::string;

#ifndef _SOLUTION_H
#define _SOLUTION_H

struct Solution
{
	double **S;//a set of solutions
	double **IniS;//initial solutions, more than S
	double **UsedIniS;//record the used Initial solutions, for sensitivity testing
	double **OriginalIniS;//record the first time generated initial solutions
	double *GlobalOptima;//record the position of the global optima
	
	//double seed;//random seed to initiate solutions
	//int NumSeed;//number of seed tested
	double *RangeofNeighborhood;//define how large is the range of neighborhood
	//double *MaxRangeNeighbor;
	int NumRange;//define how many neighborhood ranges will be tested
	int NumMove;//define the number of initial solutions that will be moved into the neighborhood

	double *AvgConvergence;//convergence recording
	int NumConvergence;//the number of convergence point to be recorded
	int NumConvergenceStep;//the step of each convergence record, say, 10 FE record once
	double BestEver;
	double *LoopResults;
	string NeighborMethod;//to select the neighbor moving method used, gaussian or square
	double MinNeighborRange;//min and max range of neighborhood
	double MaxNeighborRange;//
	int NeighborCount;// number of neighbor ranges to be tested

	double **MeanperMoveNeighbor;//to record the mean for each moved and each neighbor setting
	double **StdperMoveNeighbor;//to record the std for each moved and each neighbor setting

	int EANP;//fixed size initial population size for EA
	int Pdouble;
	int IniNP;//number of solutions initialized by techniques, can be larger than P
	double *Min;
	double *Max;//min and max value of the variable
	int D;//number of variables in each solution, dimension
	double *Y;//objective function value of the solution
	double *IniY;//objective function value of the initial solutions
	int Func_num;
	
	int Tech_num;//1 for chaos, 2 for opposition
	double TargetBest;//theoritically optima point
	double Error;//if the difference between the current best and TargetBest is smaller than this Error, stop and record FE
	//int *FEused;//to record how many FE used to reach the target range
	bool AdaptiveFlag;//to turn on/off the adaptive parameter updating scheme for each algorithm
	char *dir_file;//the file recording the sobol directions

	int Algorithm;//name of the optimization algorithm used
	string ConstraintHandle;//name of the method used to handle the constraints violation

	int MaxFE;//maximum FE times
	int Looptime;//loop time

	int CECorBBOB;//to decide whether CEC or BBOB benchmark functions to use
	int trialID;//to decide the rseed used for BBOB09
	int CRODecomS;
};

class SolutionOperator
{
public:
	SolutionOperator();
	virtual~SolutionOperator();

	void UpdateFunction(int, int, int);
	void UpdateY();
	void Release(int);
};

#endif