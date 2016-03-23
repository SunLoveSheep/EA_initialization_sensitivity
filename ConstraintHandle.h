//Functions to handle the constraint violations

#ifndef _CONSTRAINTHANDLE_H
#define _CONSTRAINTHANDLE_H

class ConstraintHandle
{
public:
	ConstraintHandle();
	virtual~ConstraintHandle();

	void Equaltobound(double *X);//when violating, equal the number to the boundary
	void Bounceback(double *X);//when violating, bounce the number back by the boundary

	void Handleconstaint(double *X);
};

#endif