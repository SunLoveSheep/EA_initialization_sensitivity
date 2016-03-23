#ifndef _cec14_H
#define _cec14_H

class cec14{
public:
	cec14();
	virtual ~cec14();

	void FileIO(int N, int func_num);
	void Release();

	double cec14_problems(const double * x, int dim, int no);
};

#endif