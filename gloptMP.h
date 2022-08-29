#pragma once

#include "hilbert.h"

class GloptMP
{
public:
	int ndim; // Number of dimensions
	int nbits; // Number of bits for each dimensions
	
	int nvcur, nvmax, nproc, niter, iiter, stop;
	double rel, eps_equal, eps_stop, x_min, f_min;
	
	HilbertCurve *hc;
	std::vector <std::pair <double,double>> TrialPoints;
	
	GloptMP (int _ndim, int _nbits, int _nproc, int _niter, double _rel, double _eps_equal);
	~GloptMP ();
	
	void getNextPoints (double *coords); // returns an array of coords [ndim * nproc]
	void pushNextPoint (double x, double func);
};
