#pragma once

#include <cstdio>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

typedef unsigned int coord_t; // char, short, int for up to 8,16,32 bits per word

class HilbertCurve
{
public:
	int ndim; // Number of dimensions
	int nbits; // Number of bits for each dimensions
	
	HilbertCurve (int _ndim, int _nbits) : ndim (_ndim), nbits (_nbits) {};
	
	void Trans2Axes (coord_t* X); // Convert Hilbert transposed coordinates to Cartesian Z^n coordinates
	void Axes2Trans (coord_t* X); // Convert Cartesian Z^n coordinates to Hilbert transposed coordinates 
	unsigned long long Trans2Long (coord_t* X); // Convert Hilbert transposed coordinates to a long integer number
	void Long2Trans (coord_t* X, unsigned long long code); // Convert a long integer number to Hilbert transposed coordinates
	void Double2DoubleAxes (double* X, double code, int scalein = 1, int scaleout = 0); // Convert a long integer number to Cartesian R^n coordinates
	void DrawHilbert (std::vector <std::pair <double,double>> TrialPoints, int nv); // Draw a 2D Hilbert curve with a set of points
};
