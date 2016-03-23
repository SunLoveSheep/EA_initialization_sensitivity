#include "stdafx.h"
#include "CRO.h"
#include "ConstraintHandle.h"
/*#include "solution.h"
#include "Initialization.h"
#include "ConstraintHandle.h"
#include "FEoutput.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include "cec14_eotp.h"
#include "cec13.h"
#include "cec14.h"
#include "bbob09functions.h"*/


using std::string;
using namespace std;

CROparameter CROp;

CRO::CRO()
{

}

CRO::~CRO()
{
	//delete []CROp.NumHit;
	//delete []CROp.MinHit;
}

void CRO::CROpInitialization()
{
	extern Solution solution;
	
	//CROp.InitialKE=10;
	CROp.KElossrate=0.7;
	CROp.MCollisionrate=0.2;//t>this value, single molecule, t<this value, multi
	CROp.SynReduce=1;//=1, no reduce, >1 reduce
	//CROp.Beta=CROp.InitialKE/20;
	CROp.Alpha=(solution.MaxFE/10)/solution.EANP;//more P, less alpha.
	CROp.NumHit=new double[solution.Pdouble];
	CROp.MinHit=new double[solution.Pdouble];
	for (int p=0;p<solution.Pdouble;p++)
	{
		CROp.NumHit[p]=0;
		CROp.MinHit[p]=0;
	}
}

void CRO::Gaussneighbor(double *In, double *Out)
{
	extern Solution solution;
	Initialization initialization;
	ConstraintHandle constrainthandle;
	
	for (int i=0;i<solution.D;i++)
	{
		Out[i]=initialization.GaussRandomGen(In[i],CROsigma);
	}

	constrainthandle.Handleconstaint(Out);
}

void CRO::Synthesisneighbor(double *Out, double *X1, double *X2)
{
	extern Solution solution;
	Initialization initialization;
	ConstraintHandle constrainthandle;
	
	double p=0;

	for (int i=0;i<solution.D;i++)
	{
		p=initialization.RandomGen(0,1);
		if (p<0.5)
		{
			Out[i]=X1[i];
		}
		else
		{
			Out[i]=X2[i];
		}
	}

	constrainthandle.Handleconstaint(Out);
}

void CRO::CROonwall(double *X, int *M)
{
	extern Solution solution;
	Initialization initialization;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	double *temp=new double[solution.D];
	double PEresult;
	double Energyjudge=0;//to judge whether the reation will happen
	double q;//a random number between 0 and 1

	//temp=Gaussneighbor(X);
	Gaussneighbor(X,temp);
	
	if (solution.CECorBBOB==0)
		PEresult=CEC14.cec14_eotp_problems(temp,solution.D,solution.Func_num);
	else if (solution.CECorBBOB==1)
		PEresult=bbob09.FunctionCalculation(temp,solution.Func_num);
	else if (solution.CECorBBOB==2)
		PEresult=CEC13.cec13_problems(temp,solution.D,solution.Func_num);
	else if (solution.CECorBBOB==3)
		PEresult=CEC14_normal.cec14_problems(temp,solution.D,solution.Func_num);

	FE++;//updating function evaluation time
	CROp.NumHit[M[0]]++;
	Energyjudge=PE[M[0]]+KE[M[0]]-PEresult;

	if ((PE[M[0]]<0)&&(KE[M[0]]<0)&&(PEresult<0))
		Energyjudge=Energyjudge*(-1);

	if (Energyjudge>=0)//onwall happens
	{
		q=initialization.RandomGen(CROp.KElossrate,1);
		KE[M[0]]=Energyjudge*q;
		PE[M[0]]=PEresult;
		solution.Y[M[0]]=PEresult;
		buffer=buffer+Energyjudge*(1-q);
		
		//update the solution information
		for (int i=0;i<solution.D;i++)
		{
			solution.S[M[0]][i]=temp[i];
		}

		CROp.MinHit[M[0]]=CROp.NumHit[M[0]];
		Onwallcounter++;
	}

	delete []temp;
}

void CRO::CROdecomposition(double *X, int *M, bool *J)
{
	extern Solution solution;
	Initialization initialization;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	if (CROp.NumHit[M[0]]-CROp.MinHit[M[0]]<=CROp.Alpha)
	{
		return;
	}
	
	double *temp1=new double[solution.D];
	double *temp2=new double[solution.D];
	double KEresult1,KEresult2,PEresult1,PEresult2;
	double Energyjudge=0;//to decide whether the decom will happen
	double k=0;//a random number between 0 and 1 to update energies

	Gaussneighbor(X,temp1);
	Gaussneighbor(X,temp2);

	if (solution.CECorBBOB==0)
	{
		PEresult1=CEC14.cec14_eotp_problems(temp1,solution.D,solution.Func_num);
		PEresult2=CEC14.cec14_eotp_problems(temp2,solution.D,solution.Func_num);
	}
	else if (solution.CECorBBOB==1)
	{
		PEresult1=bbob09.FunctionCalculation(temp1,solution.Func_num);
		PEresult2=bbob09.FunctionCalculation(temp2,solution.Func_num);
	}
	else if (solution.CECorBBOB==2)
	{
		PEresult1=CEC13.cec13_problems(temp1,solution.D,solution.Func_num);
		PEresult2=CEC13.cec13_problems(temp2,solution.D,solution.Func_num);
	}
	else if (solution.CECorBBOB==3)
	{
		PEresult1=CEC14_normal.cec14_problems(temp1,solution.D,solution.Func_num);
		PEresult2=CEC14_normal.cec14_problems(temp2,solution.D,solution.Func_num);
	}
	FE=FE+2;

	//Energyjudge=cro.PE[M[0]]+cro.KE[M[0]]-PEresult1-PEresult2;
	Energyjudge=PE[M[0]]+KE[M[0]]-PEresult1-PEresult2;

	if ((PE[M[0]]<0)&&(KE[M[0]]<0)&&(PEresult1<0)&&(PEresult2<0))
		Energyjudge=Energyjudge*(-1);

	//cout<<Energyjudge<<" "<<buffer<<endl;;
	if (Energyjudge>=0)
	{
		J[0]=1;//update the judger, the decom happens
		k=initialization.RandomGen(0,1);
		KEresult1=KE[M[0]]*k;
		KEresult2=KE[M[0]]*(1-k);

		solution.EANP++;
		//to prevent overheading
		if (solution.EANP>=solution.Pdouble)
		{
			cout<<"warning!! Population size has been doubled!"<<endl;
			FE=FE-2;
			return;
		}
		//update solution space
		for (int i=0;i<solution.D;i++)
		{
			solution.S[M[0]][i]=temp1[i];
			solution.S[solution.EANP-1][i]=temp2[i];
		}
		KE[M[0]]=KEresult1;
		KE[solution.EANP-1]=KEresult2;
		PE[M[0]]=PEresult1;
		PE[solution.EANP-1]=PEresult2;
		solution.Y[M[0]]=PEresult1;
		solution.Y[solution.EANP-1]=PEresult2;

		CROp.NumHit[M[0]]=0;
		CROp.NumHit[solution.EANP-1]=0;
		CROp.MinHit[M[0]]=0;
		CROp.MinHit[solution.EANP-1]=0;
		Decomcounter++;
	}
	else if (Energyjudge+buffer>=0)
	{
		J[0]=1;//update the judger, the decom happens
		double m1,m2,m3,m4;
		m1=initialization.RandomGen(0,1);
		m2=initialization.RandomGen(0,1);
		m3=initialization.RandomGen(0,1);
		m4=initialization.RandomGen(0,1);
		KEresult1=(Energyjudge+buffer)*m1*m2;
		KEresult2=(Energyjudge+buffer-KEresult1)*m3*m4;
		buffer=buffer+Energyjudge-KEresult1-KEresult2;

		solution.EANP++;
		//to prevent overheading
		if (solution.EANP>=solution.Pdouble)
		{
			cout<<"warning!! Population size has been doubled!";
			FE=FE-2;
			return;
		}

		for (int i=0;i<solution.D;i++)
		{
			solution.S[M[0]][i]=temp1[i];
			solution.S[solution.EANP-1][i]=temp2[i];
		}
		KE[M[0]]=KEresult1;
		KE[solution.EANP-1]=KEresult2;
		PE[M[0]]=PEresult1;
		PE[solution.EANP-1]=PEresult2;
		solution.Y[M[0]]=PEresult1;
		solution.Y[solution.EANP-1]=PEresult2;

		CROp.NumHit[M[0]]=0;
		CROp.NumHit[solution.EANP-1]=0;
		CROp.MinHit[M[0]]=0;
		CROp.MinHit[solution.EANP-1]=0;
		Decomcounter++;
	}
	else
	{
		J[0]=0;
	}

	delete []temp1;
	delete []temp2;
}

void CRO::CROintermolecularcollision(double *X1, double *X2, int *M)
{
	extern Solution solution;
	Initialization initialization;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	double *temp1=new double[solution.D];
	double *temp2=new double[solution.D];
	double KEresult1, KEresult2, PEresult1, PEresult2;
	double Energyjudge=0;
	double p;//random number between 0,1

	Gaussneighbor(X1,temp1);
	Gaussneighbor(X2,temp1);
		
	if (solution.CECorBBOB==0)
	{
		PEresult1=CEC14.cec14_eotp_problems(temp1,solution.D,solution.Func_num);
		PEresult2=CEC14.cec14_eotp_problems(temp2,solution.D,solution.Func_num);
	}
	else if (solution.CECorBBOB==1)
	{
		PEresult1=bbob09.FunctionCalculation(temp1,solution.Func_num);
		PEresult2=bbob09.FunctionCalculation(temp2,solution.Func_num);
	}
	else if (solution.CECorBBOB==2)
	{
		PEresult1=CEC13.cec13_problems(temp1,solution.D,solution.Func_num);
		PEresult2=CEC13.cec13_problems(temp2,solution.D,solution.Func_num);
	}
	else if (solution.CECorBBOB==3)
	{
		PEresult1=CEC14_normal.cec14_problems(temp1,solution.D,solution.Func_num);
		PEresult2=CEC14_normal.cec14_problems(temp2,solution.D,solution.Func_num);
	}

	FE=FE+2;
	CROp.NumHit[M[0]]++;
	CROp.NumHit[M[1]]++;
	Energyjudge=(PE[M[0]]+KE[M[0]]+PE[M[1]]+KE[M[1]])-(PEresult1+PEresult2);

	if ((PE[M[0]]<0)&&(KE[M[0]]<0)&&(PEresult1<0)&&(PEresult2<0)&&(PE[M[1]]<0)&&(KE[M[1]]<0))
		Energyjudge=Energyjudge*(-1);

	if (Energyjudge>=0)
	{
		p=initialization.RandomGen(0,1);
		KEresult1=Energyjudge*p;
		KEresult2=Energyjudge*(1-p);

		for (int i=0;i<solution.D;i++)
		{
			solution.S[M[0]][i]=temp1[i];
			solution.S[M[1]][i]=temp2[i];
		}
		KE[M[0]]=KEresult1;
		KE[M[1]]=KEresult2;
		PE[M[0]]=PEresult1;
		PE[M[1]]=PEresult2;
		//solution.Y[M[0]]=PEresult1;
		//solution.Y[M[1]]=PEresult2;

		CROp.MinHit[M[0]]=CROp.NumHit[M[0]];
		CROp.MinHit[M[1]]=CROp.NumHit[M[1]];
		Intercounter++;
	}

	delete []temp1;
	delete []temp2;
}

void CRO::CROsynthesis(double *X1, double *X2, int *M, bool *J)
{
	extern Solution solution;
	Initialization initialization;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	if ((KE[M[0]]>=CROp.Beta)||(KE[M[1]]>=CROp.Beta))
	{
		return;
	}
	
	double *temp=new double[solution.D];
	double KEresult,PEresult;
	double Energyjudge=0;

	//temp=Synthesisneighbor(X1,X2);
	Synthesisneighbor(temp,X1,X2);
	
	if (solution.CECorBBOB==0)
	{
		PEresult=CEC14.cec14_eotp_problems(temp,solution.D,solution.Func_num);
	}
	else if (solution.CECorBBOB==1)
	{
		PEresult=bbob09.FunctionCalculation(temp,solution.Func_num);
	}
	else if (solution.CECorBBOB==2)
	{
		PEresult=CEC13.cec13_problems(temp,solution.D,solution.Func_num);
	}
	else if (solution.CECorBBOB==3)
	{
		PEresult=CEC14_normal.cec14_problems(temp,solution.D,solution.Func_num);
	}

	FE++;
	Energyjudge=PE[M[0]]+KE[M[0]]+PE[M[1]]+KE[M[1]]-CROp.SynReduce*PEresult;
	//Synreduce to reduce the possibility of a synthesis

	if ((PE[M[0]]<0)&&(KE[M[0]]<0)&&(PEresult<0)&&(PE[M[1]]<0)&&(KE[M[1]]<0))
		Energyjudge=Energyjudge*(-1);

	if (Energyjudge>=0)
	{
		J[0]=1;
		KEresult=Energyjudge;
		
		//copy the resultant molecule to the solution set
		for (int i=0;i<solution.D;i++)
		{
			solution.S[M[0]][i]=temp[i];
		}
		PE[M[0]]=PEresult;
		KE[M[0]]=KEresult;

		//reduce the other reactant molecule from the solution set
		for (int p=M[1]+1;p<solution.EANP;p++)
		{
			for (int n=0;n<solution.D;n++)
			{
				solution.S[p-1][n]=solution.S[p][n];
			}
			KE[p-1]=KE[p];
			PE[p-1]=PE[p];
			solution.Y[p-1]=solution.Y[p];
		}

		solution.EANP--;
		CROp.MinHit[M[0]]=0;
		CROp.NumHit[M[0]]=0;
		Syncounter++;
	}

	delete []temp;
}

void CRO::CROprocess(int loop)
{
	extern Solution solution;
	SolutionOperator SOperator;
	Initialization initialization;
	FEoutput feoutput;
	cec14_eotp CEC14;
	cec13 CEC13;
	cec14 CEC14_normal;
	BBOB09 bbob09;
	
	feoutput.Initial();
	//SOperator.UpdateY();
	double min=1000000000000000000;
	for (int i=0;i<solution.EANP;i++)
	{
		if (solution.Y[i]<min)
		{
			min=solution.Y[i];
		}
	}
	CROpInitialization();
	CROp.InitialKE=min/4;
	CROp.Beta=CROp.InitialKE/20;

	double *X1=new double[solution.D];//to temporarily store one molecule
	double *X2=new double[solution.D];//to temporarily store one molecule
	KE=new double[solution.Pdouble];
	PE=new double[solution.Pdouble];
	double t;//to check with the molecular collision rate
	int *M=new int[2];//to choose molecules from the population, M[0] for the 1st picked molecule. M[1] for the 2nd
	M[0]=M[1]=0;
	bool *J=new bool[1];//to check if a Decom or Synth has happened. J=1 means happened, J=0 means not happened
	J[0]=0;
	
	CROsigmaAdaptive=0.1*(solution.Max[0]-solution.Min[0]);
	CROsigmaNonAdaptive=0.05*(solution.Max[0]-solution.Min[0]);
	if (solution.AdaptiveFlag=1)
	{
		CROsigma=CROsigmaAdaptive;
	}
	else
	{
		CROsigma=CROsigmaNonAdaptive;
	}
	double Sigma=CROsigma;
	double FEmax=solution.MaxFE;

	//initial the molecules' KE and PE
	for (int p=0;p<solution.EANP;p++)
	{
		for (int n=0;n<solution.D;n++)
		{
			X1[n]=solution.S[p][n];
		}
		
		//KE[p]=CROp.InitialKE;
		if (solution.CECorBBOB==0)
		{
			PE[p]=CEC14.cec14_eotp_problems(X1,solution.D,solution.Func_num);
		}
		else if (solution.CECorBBOB==1)
		{
			PE[p]=bbob09.FunctionCalculation(X1,solution.Func_num);
		}
		else if (solution.CECorBBOB==2)
		{
			PE[p]=CEC13.cec13_problems(X1,solution.D,solution.Func_num);
		}
		else if (solution.CECorBBOB==3)
		{
			PE[p]=CEC14_normal.cec14_problems(X1,solution.D,solution.Func_num);
		}
		//PE[p] = solution.Y[p];
	}
	
	buffer=0;
	FE=0;
	Onwallcounter=0;
	Decomcounter=0;
	Intercounter=0;
	Syncounter=0;

	do 
	{
		if (solution.AdaptiveFlag==1)
		{
			//self adaptive of stepsize
			CROsigma=Sigma-Sigma*(double(FE)/double(FEmax));//less aggresive
			//CROsigma=CROsigma-CROsigma*(FE/FEmax);//more aggresive, to 0 at some FE, no more neighborsearch

			//parameter reset

		}
		
		t=initialization.RandomGen(0,1);
		J[0]=0;
		if ((t>CROp.MCollisionrate)||(solution.EANP==1))//single molecule reaction happen
		{
			M[0]=rand()%(solution.EANP);//randomly pick one solution from the population
			for (int i=0;i<solution.D;i++)
			{
				X1[i]=solution.S[M[0]][i];
			}

			CROdecomposition(X1,M,J);//decomposition
			if (J[0]==0)
			{
				CROonwall(X1,M);//onwall
			}
		}
		else
		{
			M[0]=rand()%(solution.EANP);
			do
			{
				M[1]=rand()%(solution.EANP);
			}while(M[1]==M[0]);//to generate a M[1] different from M[0]

			//copy the to selected molecules
			for (int i=0;i<solution.D;i++)
			{
				X1[i]=solution.S[M[0]][i];
				X2[i]=solution.S[M[1]][i];
			}

			CROsynthesis(X1,X2,M,J);
			if (J[0]==0)
			{
				CROintermolecularcollision(X1,X2,M);
			}
		}
		
		//SOperator.UpdateY();
		//feoutput.FEOprocess(FE);
		feoutput.UpdateBest(FE);
	}while(FE<solution.MaxFE);
	//end of while loop
	/*
	cout<<"Reaction Counters: "<<endl;
	cout<<"Onwall: "<<Onwallcounter<<endl;
	cout<<"Decomposition: "<<Decomcounter<<endl;
	cout<<"Intermolecular: "<<Intercounter<<endl;
	cout<<"Synthesis: "<<Syncounter<<endl;
	*/
	
	feoutput.Release();
	
	delete []CROp.MinHit;
	delete []CROp.NumHit;
	delete []KE;
	delete []PE;
	delete []X1;
	delete []X2;
	delete []M;
	delete []J;
}