#ifndef _bbob09support_H
#define _bbob09support_H

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

class BBOB09support{
public:
	BBOB09support();
	virtual~ BBOB09support();

	void gauss(double * g, int N, int seed);
	double round(double a);
	double fmin(double a, double b);
	double fmax(double a, double b);
	void unif(double* r, int N, int inseed);
	void computeXopt(double *xshift, int seed);
	double computeFopt(int _funcId);
	void monotoneTFosc(double* f);
	void computeRotation(double ** B, int seed);
};

#endif

