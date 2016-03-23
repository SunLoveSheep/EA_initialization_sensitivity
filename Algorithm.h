//Optimization algorithms set
#include<string>

#ifndef _ALGORITHM_H
#define _ALGORITHM_H

using std::string;

class Algorithm
{
public:
	Algorithm();
	virtual~Algorithm();

	void Optimization(int loop);//differential evolution
};

#endif