#include <WINDOWS.H>      
#include <stdio.h>
#include <math.h>
#include <malloc.h>
//#include <mex.h> //for matlab

#ifndef _cec13_H
#define _cec13_H

class cec13{
public:
	cec13();
	virtual ~cec13();

	void FileIO(int dim, int func_num);
	void Release();
	
	double cec13_problems(const double * x, int dim, int no);

	//double sphere_func (const double *, int , double *,double *, int); /* Sphere */
	//double ellips_func(const double *, int , double *,double *, int); /* Ellipsoidal */
	//double bent_cigar_func(const double *, int , double *,double *, int); /* Discus */
	//double discus_func(const double *, int , double *,double *, int);  /* Bent_Cigar */
	//double dif_powers_func(const double *, int , double *,double *, int);  /* Different Powers */
	//double rosenbrock_func (const double *, int , double *,double *, int); /* Rosenbrock's */
	///double schaffer_F7_func (const double *, int , double *,double *, int); /* Schwefel's F7 */
	//double ackley_func (const double *, int , double *,double *, int); /* Ackley's */
	//double rastrigin_func (const double *, int , double *,double *, int); /* Rastrigin's  */
	//double step_rastrigin_func (const double *, int , double *,double *, int); /* Noncontinuous Rastrigin's  */
	//double weierstrass_func (const double *, int , double *,double *, int); /* Weierstrass's  */
	//double griewank_func (const double *, int , double *,double *, int); /* Griewank's  */
	//double schwefel_func (const double *, int , double *,double *, int); /* Schwefel's */
	//double katsuura_func (const double *, int , double *,double *, int); /* Katsuura */
	//double bi_rastrigin_func (const double *, int , double *,double *, int); /* Lunacek Bi_rastrigin Function */
	//double grie_rosen_func (const double *, int , double *,double *, int); /* Griewank-Rosenbrock  */
	//double escaffer6_func (const double *, int , double *,double *, int); /* Expanded Scaffer¡¯s F6  */
	//double cf01 (const double *, int , double *,double *, int); /* Composition Function 1 */
	//double cf02 (const double *, int , double *,double *, int); /* Composition Function 2 */
	//double cf03 (const double *, int , double *,double *, int); /* Composition Function 3 */
	//double cf04 (const double *, int , double *,double *, int); /* Composition Function 4 */
	//double cf05 (const double *, int , double *,double *, int); /* Composition Function 5 */
	//double cf06 (const double *, int , double *,double *, int); /* Composition Function 6 */
	//double cf07 (const double *, int , double *,double *, int); /* Composition Function 7 */
	//double cf08 (const double *, int , double *,double *, int); /* Composition Function 8 */

	//void shiftfunc (const double*,double*,int,double*);
	//void rotatefunc (const double*,double*,int, double*);
	//void asyfunc (const double *, double *x, int, double);
	//void oszfunc (const double *, double *, int);
	//void cf_cal(const double *, double *, int, double *,double *,double *,double *,int);
};

#endif