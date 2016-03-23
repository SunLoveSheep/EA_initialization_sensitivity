#include "stdafx.h"
#include "ConstraintHandle.h"
//#include "solution.h"

ConstraintHandle::ConstraintHandle()
{
	
}

ConstraintHandle::~ConstraintHandle()
{

}

void ConstraintHandle::Equaltobound(double *X)
{
	extern Solution solution;

	for (int i=0;i<solution.D;i++)
	{
		if(X[i]<solution.Min[i])
		{
			X[i]=solution.Min[i];
		}
		else if(X[i]>solution.Max[i])
		{
			X[i]=solution.Max[i];
		}
	}

}

void ConstraintHandle::Bounceback(double *X)
{
	extern Solution solution;

	for (int i=0;i<solution.D;i++)
	{
		if(X[i]<solution.Min[i])
		{
			X[i]=2*solution.Min[i]-X[i];
		}
		else if(X[i]>solution.Max[i])
		{
			X[i]=2*solution.Max[i]-X[i];
		}
	}
}

void ConstraintHandle::Handleconstaint(double *X)
{
	extern Solution solution;

	if (solution.ConstraintHandle=="Equaltobound")
	{
		Equaltobound(X);
	}
	
	if (solution.ConstraintHandle=="Bounceback")
	{
		Bounceback(X);
	}

}