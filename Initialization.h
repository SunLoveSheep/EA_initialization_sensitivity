//Initialization functions
#include <string>

#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H

using std::string;

class Initialization
{
public:
	Initialization();
	virtual~Initialization();

	void Initial();//initial the solutions
	double RandomGen(double min, double max);//generate a random number between min and max
	double GaussRandomGen(double miu, double sigma);//generate a gaussian random number

	void MoveIniSolution(int);//move initial solutions into the neighborhood

	void Random();//chaos technique
	void Chaos();//chaos technique
	void Opposition();//opposition technique
	void QuasiOpposition();//quasi-opposition technique
	void QuasiInterpolation();//quasi interpolation technique
	void Sobol();//sobol technique
};

#endif