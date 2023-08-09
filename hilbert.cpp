#include "hilbert.h"

//+++++++++++++++++++++++++++ PUBLIC-DOMAIN SOFTWARE ++++++++++++++++++++++++++
// Functions: trans2Axes axes2Trans
// Purpose: Transform in-place between Hilbert transpose and geometrical axes
// Example: b=5 bits for each of n=3 coordinates.
// 15-bit Hilbert integer = A B C D E F G H I J K L M N O is stored
// as its Transpose
//        X[0] = A D G J M             X[2]|
//        X[1] = B E H K N <------->       | /X[1]
//        X[2] = C F I L O            axes |/
//               high  low                 0------ X[0]
// Axes are stored conventially as b-bit integers.
// Author: John Skilling 20 Apr 2001 to 11 Oct 2003
//-----------------------------------------------------------------------------

// Convert Hilbert transposed coordinates to Cartesian Z^n coordinates
void HilbertCurve::Trans2Axes (coord_t* X)
{
    coord_t N = 2 << (nbits - 1), P, Q, t;

    // Gray decode by H ^ (H/2)
    t = X[ndim - 1] >> 1;
    // Corrected error in Skilling's paper on the following line. The appendix had i >= 0 leading to negative array index.
    for (int i = ndim - 1; i > 0; i--) X[i] ^= X[i - 1];
    X[0] ^= t;

    // Undo excess work
    for (Q = 2; Q != N; Q <<= 1) {
        P = Q - 1;
        for (int i = ndim - 1; i >= 0; i--)
            if (X[i] & Q) // Invert
                X[0] ^= P;
            else { // Exchange
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
    }
}

// Convert Cartesian Z^n coordinates to Hilbert transposed coordinates 
void HilbertCurve::Axes2Trans (coord_t* X)
{
    coord_t M = 1 << (nbits - 1), P, Q, t;

    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1) {
        P = Q - 1;
        for (int i = 0; i < ndim; i++)
            if (X[i] & Q) // Invert
                X[0] ^= P;
            else { // Exchange
                t = (X[0] ^ X[i]) & P;
                X[0] ^= t;
                X[i] ^= t;
            }
    }

    // Gray encode
    for (int i = 1; i < ndim; i++) X[i] ^= X[i - 1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if (X[ndim - 1] & Q) t ^= Q - 1;
    for (int i = 0; i < ndim; i++) X[i] ^= t;
}

// Convert Hilbert transposed coordinates to a long integer number
unsigned long long HilbertCurve::Trans2Long (coord_t* X)
{
    const int nbits2 = ndim * nbits;
	unsigned int code[ndim];
	
    for (int i = 0, andbit = 1; i < nbits2; i += ndim, andbit <<= 1) {
		for (int j = 0; j < ndim; j++)
			code[j] |= (unsigned int)(X[0] & andbit) << i;
    }

    unsigned retcode = code[0] << (ndim - 1);
		for (int j = 1; j < ndim; j++)
			retcode |= (code[j] << (ndim - j - 1));
    return retcode;
}

// Convert a long integer number to Hilbert transposed coordinates
void HilbertCurve::Long2Trans (coord_t* X, unsigned long long code)
{
	for (int j = 0; j < ndim; j++)
		X[j] = 0;

    for (unsigned int i = 0; i <= nbits; ++i) {
		unsigned int selector = 1;
		unsigned int shift_selector = ndim * i;
		unsigned int shiftback = (ndim - 1) * i;
		for (int j = 0; j < ndim; j++)
			X[ndim - j - 1] |= (code & (selector << (shift_selector + j))) >> (shiftback + j);
    }
}

// Convert a double number to Cartesian R^n coordinates,
// scalin = 1 (default) means that input number has to be scaled from [0,1] to [0,2^(ndim*nbits)-1],
// scalout = 1 (default) means that output X[ndim] has to be scaled from [0,2^nbits-1]^ndim to [0,1]^ndim
void HilbertCurve::Double2DoubleAxes (double* X, double code, int scalein, int scaleout)
{
	coord_t ux1[ndim], ux2[ndim];
	double step;
	unsigned long long lcode;
	
	if (scalein == 1)
		code *= ((1 << ndim * nbits) - 1);
	
	if (code < 0)
	{
		std::cout << "Correcting negative point!" << std::endl;
		code *= -1;
	}
	
	lcode = (unsigned long long) code;
	Long2Trans (ux1, lcode);
	Long2Trans (ux2, lcode + 1);
	
	Trans2Axes (ux1);
	Trans2Axes (ux2);
	
	step = code - (double) lcode;
	for (int i = 0; i < ndim; ++i)
		X[i] = ux1[i] * (1.0 - step) + ux2[i] * step; 
	
	if (scaleout == 1)
		for (int i = 0; i < ndim; ++i)
			X[i] = (X[i] + 0.5) / (1 << nbits); 
}

void HilbertCurve::DrawHilbert (std::vector <std::pair <double,double>> TrialPoints, int nv)
{
	if (ndim == 2)
	{
		coord_t X[ndim];
		double Y[ndim];
		std::ofstream fs;
		
		fs.open ("curve.ps");
	
		fs << "%!PS-Adobe-2.0 EPSF-2.0" << std::endl;
		fs << "%%BoundingBox: 0 0 " << (1 << nbits) + 1 << " " << (1 << nbits) + 1 << std::endl;
		fs << "%%EndComments" << std::endl;
		
		fs << "/m{moveto} def" << std::endl;
		fs << "/l{lineto} def" << std::endl;

		fs << "/slw{setlinewidth}def" << std::endl;

		fs << "0 slw" << std::endl;
		fs << "1 1 translate" << std::endl;
		
		fs << "newpath" << std::endl;
		fs << 0 << " " << 0 << " m " << std::endl;
		for (unsigned long long code = 1; code < (1 << ndim * nbits); ++code)
		{
			Long2Trans (X, code);
			Trans2Axes (X);
			fs << X[0] << " " << X[1] << " l " << std::endl;
		}
		fs << "stroke" << std::endl;
		
		for (int i = 0; i < nv; ++i)
		{
			Double2DoubleAxes (Y, TrialPoints[i].first, 1, 0);
			fs << "newpath" << std::endl;
			fs << Y[0] - 0.3 << " " << Y[1] - 0.3 << " m " << std::endl;
			fs << Y[0] + 0.3 << " " << Y[1] + 0.3 << " l " << std::endl;
			fs << "stroke" << std::endl;
			fs << "newpath" << std::endl;
			fs << Y[0] - 0.3 << " " << Y[1] + 0.3 << " m " << std::endl;
			fs << Y[0] + 0.3 << " " << Y[1] - 0.3 << " l " << std::endl;
			fs << "stroke" << std::endl;
		}
		
		fs << "showpage" << std::endl;
		
		fs.close ();
	}
	/*
	else
	{
		coord_t X[ndim];
		for (unsigned int code = 1; code < (1 << ndim * nbits); ++code)
		{
			Long2Trans (X, code);
			Trans2Axes (X);
			std::cout << "{" << X[0];
			for (int i = 1; i < ndim; i++)
				std::cout << ", " << X[i];
			std::cout << "}" << std::endl;
		}
	}
	*/
}
