//Record the best ever results and output

#ifndef _FEOUTPUT_H
#define _FEOUTPUT_H

struct FEvariable
{
};

class FEoutput
{
public:
	FEoutput();
	virtual~FEoutput();

	//variables for fe loop
	double *Bestever;//to store the best solution ever found
	int BesteverFE;
	int ConvergeCounter;//How mant points we want in convergence record
	double *Convergence;//to store the converge data
	double *ConvergenceFirst10;//to store the converge data of the beginning 10% FEs
	
	//functions for fe loop
	void UpdateBest(int FE);//update the best ever solution
	//void RecordConvergence(int FE);//record the convergence data
	void Output(int FE);//output result

	//void FEOprocess(int FE);//whole update, record and output process
	
	//variables for algorithm loop
	double Loopbest,Loopworst,Loopmean,Loopstd;//best, worst, mean and std for loop results
	//double *Loopresults;

	//functions for algorithm loop
	void Initial();
	void Release();
	void Loopsort();//get best and worst
	void Loopmeanstd(int, int);//calculate mean and std
	void LOutput();
	void LOutputForR();
	//void LoopOprocess();//whole loop static and output
};

#endif