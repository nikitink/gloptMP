#include "gloptMP.h"

GloptMP::GloptMP (int _ndim, int _nbits, int _nproc, int _niter, double _rel, double _eps_equal) :
	ndim (_ndim), nbits (_nbits), nproc (_nproc), niter (_niter), rel (_rel), eps_equal (_eps_equal)
{ 
	hc = new HilbertCurve (ndim, nbits);
}
	
GloptMP::~GloptMP ()
{ 
	delete hc;
}

void GloptMP::getNextPoints (double *coords) // returns an array of coords [ndim * nproc]
{
}

void GloptMP::pushNextPoint (double x, double func)
{
	TrialPoints.push_back (std::make_pair (x, func));
}
