//CRO parameters and operators
#ifndef _CRO_H
#define _CRO_H

struct CROparameter
{
	//parameters of DE
	double InitialKE;
	double KElossrate;//ke loss rate
	double MCollisionrate;//molecular collision rate
	double SynReduce;//a coefficent added in syn energyjudge to reduce possibility of synthesis
	double Beta;//synthesis criterion threshold
	double Alpha;//decomposition criterion threshold
	double *NumHit,*MinHit;//to record number of hit and min number of hit
};

class CRO
{
public:
	CRO();
	virtual~CRO();
	
	double buffer;//energy buffer
	double *KE,*PE;//to store the KE and PE for each molecule
	int FE;//functional evaluation counter
	double CROsigma;//stepsize
	double CROsigmaAdaptive;
	double CROsigmaNonAdaptive;

	int Onwallcounter;
	int Decomcounter;
	int Intercounter;
	int Syncounter;

	void CROpInitialization();
	void Gaussneighbor(double *In, double *Out);//use gaussian random to find a neighborhood
	void Synthesisneighbor(double *Out, double *X1, double *X2);//synthesis operator, combining two molecules

	void CROonwall(double *X, int *M);
	void CROdecomposition(double *X, int *M, bool *J);
	void CROintermolecularcollision(double *X1, double *X2, int *M);
	void CROsynthesis(double *X1, double *X2, int *M, bool *J);

	void CROprocess(int loop);//whole process of DE
};

#endif