//DE parameters and operators
#ifndef _DE_H
#define _DE_H

struct DEparameter
{
	//parameters of DE
	double F;//x1+F(x2-x3)
	double CR;//crossover rate
};

class DE
{
public:
	DE();
	virtual~DE();
	
	//X is the selected solution candidate, *R only has one element ==R[0], to record 
	//which candidate solution is selected
	void DEmutation(double *X, int *R);//mutation operator XR=Xr1+F(Xr2-Xr3)
	void DEcrossover(double *X, int *R);//crossover operator, <CR,X=XR, >CR, X=X
	void DEselection(double *X, int *R);//selection operator

	//int FE;

	void DEprocess(int loop);//whole process of DE
};

#endif