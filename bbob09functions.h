#ifndef _bbob09_H
#define _bbob09_H

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

struct twoDoubles {
double Ftrue;
double Fval;
};

typedef struct twoDoubles TwoDoubles;

/* and now the type of the benchmark functions themselves */
typedef struct twoDoubles (*bbobFunction)(double *);

struct BBOB_variable
{
	double Fopt;
	double * tmpvect;
	double *xshift;
	double * tmx;
	double ** rotation;
	double ** rot2;
	double ** linearTF;
	double * peaks21;
	double * peaks22;
	int * rperm21;
	int * rperm22;
	double ** Xlocal21;
	double ** Xlocal22;
	double ** arrScales21;
	double ** arrScales22;
	double aK[12];
    double bK[12];
    double F0;
	double * peakvalues21;
	double * peakvalues22;
};

class BBOB09{
public:
	BBOB09();
	virtual ~BBOB09();

	void BBOBparameterInitial();
	void Release();
	double FunctionCalculation(double *x, int FunNum);
	void CalculateParameter(int FunNum);

	/*double f1(double *x);
	double f2(double *x);
	double f3(double *x);
	double f4(double *x);
	double f5(double *x);
	double f6(double *x);
	double f7(double *x);
	double f8(double *x);
	double f9(double *x);
	double f10(double *x);
	double f11(double *x);
	double f12(double *x);
	double f13(double *x);
	double f14(double *x);
	double f15(double *x);
	double f16(double *x);
	double f17(double *x);
	double f18(double *x);
	double f19(double *x);
	double f20(double *x);
	double f21(double *x);
	double f22(double *x);
	double f23(double *x);
	double f24(double *x);*/
};

#endif
