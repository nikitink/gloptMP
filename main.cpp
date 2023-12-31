#include "hilbert.h"
#include "gmshfrac.h"
#include "gloptMP.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include <omp.h>
#include <unistd.h>
#include <immintrin.h>

#ifdef INMOST_3PHASE
#include "/home/nikitink/Work/filter-project/INMOST/Source/Headers/inmost.h"
#include "/home/nikitink/Work/filter-project/threephase/multiphase.h"
#endif

#ifdef VIRTUAL_SUTURER
#include "oneleaf.h"
#endif

std::ofstream STATS;

// #define DEBUG_PRINT

bool CmpBySecond (std::pair<double, int> p1, std::pair<double, int> p2) { 
	return p1.second < p2.second; 
} 

template <typename T>
void InsertMultiple (std::vector<T> & TrialPoints, std::vector<T> & NewPoints, std::vector<std::pair<double,int>> & Indices) {
	auto size = TrialPoints.size() + NewPoints.size();
	auto tmp = std::vector<T>();
	tmp.reserve(size);

	auto v = std::size_t{0}; // index into TrialPoints
	auto p = std::size_t{0}; // index into Indices / NewPoints

	for (auto i = std::size_t{0}; i != size; ) { // index into y
		if (p < Indices.size () && v == Indices[p].second) {
			tmp.push_back(NewPoints[p]);
			++p;
			++i;
		}
		else {
			tmp.push_back(TrialPoints[v]);
			++v;
			++i;
		}
	}
	TrialPoints = tmp;
}

double F1 (double *X, int ndim) {
	if (ndim != 2)
		return 0;
	
	double bbox[4] = {0, 5, 0, 7}, Y[ndim];
	for (int i = 0; i < ndim; ++i)
		Y[i] = bbox[2 * i] + X[i] * (bbox[2 * i + 1] - bbox[2 * i]);
	
	return cos (Y[0]) * sin (Y[1]);
}

// Rosenbrock function
double F2 (double *X, int ndim) {
	if (ndim < 2)
		return 0;
	double bbox[2 * ndim], Y[ndim], sum = 0.0;
	for (int i = 0; i < ndim; ++i)
		bbox[2 * i] = -25, bbox[2 * i + 1] = 25;
	for (int i = 0; i < ndim; ++i)
		Y[i] = bbox[2 * i] + X[i] * (bbox[2 * i + 1] - bbox[2 * i]);
	
	for (int i = 0; i < ndim - 1; ++i)
		sum += (1 - X[i]) * (1 - X[i]) + 100.0 * (X[i + 1] - X[i] * X[i]) * (X[i + 1] - X[i] * X[i]);
	
	return sum;
}

#ifdef INMOST_3PHASE
double F3 (double *X, int ndim) {
	static int numgrd = 0;
	if (ndim != 4) {
		std::cout << "This test requires dim = 4" << std::endl;
		exit (1);
	}
	
	if (sqrt ((X[0] - X[2]) * (X[0] - X[2]) + (X[1] - X[3]) * (X[1] - X[3])) < 0.0625)
		return 1e+5;
		
	std::stringstream fname;
	fname << "frac-" << numgrd++ << ".vtk";
	
	int ncells = gmshFracGrid (X, fname.str().c_str());
	addFracMaterial (ncells, fname.str().c_str());
	remove ("frac.vtk");
	symlink (fname.str().c_str(), "frac.vtk");
	
	mesh_params mparams;
	linear_solver_params lparams;
	nonlinear_solver_params nparams;
	time_solver_params tparams;
	Automatizator aut;
	Model model(aut);
	MultiphaseModel mmodel;
	if (read_file("params.txt", mmodel, lparams, nparams, tparams, mparams)) {
		lparams.verbosity = 1;
		nparams.verbosity = 1;
		tparams.verbosity = 1;
		mparams.LoadMesh(model);
		Tag ftag = mparams.mesh->GetTag("CellEntityIds");
		
		model.AddSubModel("multiphase", mmodel);
		model.Initialize();
		model.PrepareEntries();
		aut.EnumerateEntries();
		tparams.SolveTime(model, nparams, lparams, mparams);
		std::cout << "Production oil " << mmodel.getOilProd() << " water " << mmodel.getWaterProd() << " gas " << mmodel.getGasProd() << std::endl;
		
		mparams.mesh->Save(fname.str().c_str());
	}
	
	return 0;
}
#endif

// double D in [9, 34], start D = 19
// double H in [5, 25], start H = 11 
// double H_f in [-H/2, 6], start H_f = 2
#ifdef VIRTUAL_SUTURER
double F4 (double *X, int ndim) {
	if (ndim < 2 || ndim > 3)
		return 0;
	double bbox[2 * ndim], sum = 0.0;
	bbox[0] = 10, bbox[1] = 33;
	bbox[2] = 6, bbox[3] = 24;
	if (ndim == 3)
		bbox[4] = 2, bbox[5] = 6;
	
	SimplifiedTemplateGeometry g;
	g.D = bbox[0] + X[0] * (bbox[1] - bbox[0]), 
	g.H = bbox[2] + X[1] * (bbox[3] - bbox[2]), 
	g.H_f = 2.0;
	if (ndim == 3)
		g.H_f = bbox[4] + X[2] * (bbox[5] - bbox[4]);
	char *argv[] = {""};
	
	auto res = performOneLeafSimulation(g, TestConstantParameters(), 1, argv);	
	EvaluatedValues EVres = *reinterpret_cast<EvaluatedValues*>(&res);
	std::cout << EVres << std::endl;
	STATS << EVres << std::endl;
	
	// Functional: [10, 6]  -> 4e+6
	//             [19, 11] -> 3e+3
	//             [33, 24] -> 6e+6
	sum = (EVres.form_closed_valve        ? 1.0 : (exp (EVres.closure_degree))) *
		  (EVres.central_coaptation > 1.5 ? (1.0 + EVres.central_coaptation - 1.5) : exp (1.5 - EVres.central_coaptation)) *
		  (EVres.effective_height < 9.0   ? exp (9.0 - EVres.effective_height) : 
										    (EVres.effective_height > 11.0 ? exp (EVres.effective_height - 11.0) : 1.0)) *
		  (EVres.length_coaptation > 6.0  ? (1.0 + EVres.length_coaptation - 6.0) : exp (6.0 - EVres.length_coaptation)) *
		  (EVres.billowing < 1.4          ? 1.0 + 1.4 - EVres.billowing : exp (EVres.billowing - 1.4));
		  
	std::cout << "Multipliers: " <<
		  (EVres.form_closed_valve        ? 1.0 : (exp (EVres.closure_degree))) << " * " <<
		  (EVres.central_coaptation > 1.5 ? (1.0 + EVres.central_coaptation - 1.5) : exp (1.5 - EVres.central_coaptation)) << " * " <<
		  (EVres.effective_height < 9.0   ? exp (9.0 - EVres.effective_height) : 
										    (EVres.effective_height > 11.0 ? exp (EVres.effective_height - 11.0) : 1.0)) << " * " <<
		  (EVres.length_coaptation > 6.0  ? (1.0 + EVres.length_coaptation - 6.0) : exp (6.0 - EVres.length_coaptation)) << " * " <<
		  (EVres.billowing < 1.4          ? 1.0 + 1.4 - EVres.billowing : exp (EVres.billowing - 1.4)) << std::endl;
		  
	STATS << "Multipliers: " <<
		  (EVres.form_closed_valve        ? 1.0 : (exp (EVres.closure_degree))) << " * " <<
		  (EVres.central_coaptation > 1.5 ? (1.0 + EVres.central_coaptation - 1.5) : exp (1.5 - EVres.central_coaptation)) << " * " <<
		  (EVres.effective_height < 9.0   ? exp (9.0 - EVres.effective_height) : 
										    (EVres.effective_height > 11.0 ? exp (EVres.effective_height - 11.0) : 1.0)) << " * " <<
		  (EVres.length_coaptation > 6.0  ? (1.0 + EVres.length_coaptation - 6.0) : exp (6.0 - EVres.length_coaptation)) << " * " <<
		  (EVres.billowing < 1.4          ? 1.0 + 1.4 - EVres.billowing : exp (EVres.billowing - 1.4)) << std::endl;
		  
	if (isnan(sum) || sum > 1e+6)
		sum = 1e+6;
		  
	std::cout << "Input: D = " << g.D << ", H = " << g.H << ", H_f = " << g.H_f << ", Functional = " << sum << std::endl << std::endl;
	STATS << "Input: D = " << g.D << ", H = " << g.H << ", H_f = " << g.H_f << ", Functional = " << sum << std::endl << std::endl;
	
	// EVres.form_closed_valve; = true
	// EVres.closure_degree; > 0
	// EVres.central_coaptation; > 1.5
	// EVres.effective_height; in [9, 11]
	// EVres.length_coaptation; > 6
	// EVres.billowing; < 1.4
	// EVres.coaptation_area; // -> 0
	
	return sum;
}
#endif

double Func (HilbertCurve *hc, double x) {
	int ndim = hc->ndim;
	double Y[ndim];
	
	hc->Double2DoubleAxes (Y, x, 1, 1);
	
	// return F1 (Y, ndim);
	return F2 (Y, ndim);
	// return F3 (Y, ndim);
	// return F4 (Y, ndim);
}

int main(int argc, char **argv) {
	clock_t time_total = clock(), time_tmp;

	int ndim = 2, nvmax = 100, nbits = 5, nproc = 4, niter = 100;
	double rel = 1.1;
	double eps_equal = 0.5 / RAND_MAX, eps_stop = 1e-2;
	double x_min, f_min = 1e+12;
	int iiter, stop = 0;
	
	int option_index = 0, c;
	static const struct option opts[] = {
		{"help",   no_argument,       0, 'H'},
		{"ndim",   required_argument, 0, 'd'},
		{"nvmax",  required_argument, 0, 'v'},
		{"nbits",  required_argument, 0, 'b'},
		{"niter",  required_argument, 0, 'i'},
		{"nproc",  required_argument, 0, 'p'},
		{"rel",    required_argument, 0, 'r'},
		{"eps",    required_argument, 0, 'e'},
	};
	
	std::ostringstream help;
	help << "Usage: gloptMP [options]" << std::endl <<
			"-H,       --help         this message" << std::endl <<
			"-d NDIM,  --ndim=NDIM    problem dimension N (default = " << ndim << ")" << std::endl <<
			"-v NVMAX, --nvmax=NVMAX  starting number of points (default = " << nvmax << ")" << std::endl <<
			"-b NBITS, --nbits=NBITS  number of bits per dimension (default = " << nbits << ")" << std::endl <<
			"-i NITER, --niter=NITER  maximum number of iterations (default = " << niter << ")" << std::endl <<
			"-p NPROC, --nproc=NPROC  number of processors (default = " << nproc << ")" << std::endl <<
			"-r REL,   --rel=REL      reliability parameter of the method (default = " << rel << ")" << std::endl <<
			"-e EPS,   --eps=EPS      epsilon for the algorithm stop (default = " << eps_stop << ")" << std::endl;
	
	while (1) {
		c = getopt_long(argc, argv, "Hd:v:b:i:p:r:e:", opts, &option_index);
		if (c == -1)
			break;
		// here sscanf is used not to overwrite initial values in case of wrong input
		switch (c) {
			case 'H': std::cout << help.str(); exit(1);
			case 'd': sscanf(optarg, "%d", &ndim);  break;
			case 'v': sscanf(optarg, "%d", &nvmax); break;
			case 'b': sscanf(optarg, "%d", &nbits); break;
			case 'i': sscanf(optarg, "%d", &niter); break;
			case 'p': sscanf(optarg, "%d", &nproc); break;
			case 'r': sscanf(optarg, "%lf", &rel);  break;
			case 'e': sscanf(optarg, "%lf", &eps_stop);  break;
		}
	}
	STATS.open ("statistics.csv");
	
	GloptMP *opt = new GloptMP (ndim, nbits, nproc, niter, rel, eps_equal);
	
	double *xN_min = new double[ndim];
	
	// Initialize points on [0, 1]
	for (int i = 0; i < nvmax; ++i) {
		double x = (double) i / (nvmax - 1);
		opt->TrialPoints.push_back (std::make_pair (x, 0.0));
	}
	
	omp_set_num_threads (nproc);
// 	#pragma omp parallel for
	for (int i = 0; i < nvmax; ++i)
		opt->TrialPoints[i].second = Func (opt->hc, opt->TrialPoints[i].first);
	
#ifdef DEBUG_PRINT
		std::cout << "TrialPoints:" << std::endl;
		for (auto& x: opt->TrialPoints)
			std::cout << x.first << " / " << x.second << ", ";
		std::cout << std::endl;
#endif
		
	// Step 1: Renumber
	// Using simple sort() function to sort
	sort (opt->TrialPoints.begin (), opt->TrialPoints.end ());
	
	for (iiter = 0; iiter < niter; ++iiter) {
		double mu = 0, mui, M, dx, dz, R, x;
		
		// Step 2: Calculate magnitudes
		for (int i = 1; i < nvmax; ++i) {
			dx = (opt->TrialPoints[i].first - opt->TrialPoints[i - 1].first);
			if (dx > eps_equal) {
				dz = (opt->TrialPoints[i].second - opt->TrialPoints[i - 1].second);
				mui = fabs (dz) / pow (dx, 1.0 / ndim);
				if (mui > mu)
					mu = mui;
			}
		}
		M = (mu > 0) ? rel * mu : 1;
		
		std::vector <std::pair<double, int>> Chars;
		// Step 3: Compute characteristics
		for (int i = 1; i < nvmax; ++i) {
			R = 0;
			dx = pow (opt->TrialPoints[i].first - opt->TrialPoints[i - 1].first, 1.0 / ndim);
			dz = (opt->TrialPoints[i].second - opt->TrialPoints[i - 1].second);
			R = dx + dz * dz / (M * M * dx) - 2.0 * (opt->TrialPoints[i].second + opt->TrialPoints[i - 1].second) / M;
			Chars.push_back (std::make_pair (R, i));
		}
		
#ifdef DEBUG_PRINT
		std::cout << "1. Chars:" << std::endl;
		for (auto& x: Chars)
			std::cout << x.first << " / " << x.second << ", ";
		std::cout << std::endl;
#endif
		sort(Chars.begin(), Chars.end());
		
#ifdef DEBUG_PRINT
		std::cout << "2. Chars:" << std::endl;
		for (auto& x: Chars)
			std::cout << x.first << " / " << x.second << ", ";
		std::cout << std::endl;
#endif
		
		std::vector <std::pair<double, int>>::const_iterator first = Chars.begin() + Chars.size() - nproc;
		std::vector <std::pair<double, int>>::const_iterator last = Chars.end();
		std::vector <std::pair<double, int>> Indices(first, last);
		
		sort (Indices.begin(), Indices.end(), CmpBySecond);
		
#ifdef DEBUG_PRINT
		std::cout << "Indices:" << std::endl;
		for (auto& x: Indices)
			std::cout << x.first << " / " << x.second << ", ";
		std::cout << std::endl;
#endif
		
		std::vector <std::pair <double,double>> NewPoints;
		for (int i = 0; i < Indices.size(); ++i) {
			int iind = Indices[i].second;
			if (iind > 1) {
				dz = (opt->TrialPoints[iind].second - opt->TrialPoints[iind - 1].second);
				x = 0.5 * (opt->TrialPoints[iind - 1].first + opt->TrialPoints[iind].first) -
					((dz > 0) ? 1 : -1) * pow (fabs (dz) / mu, ndim) / (2.0 * rel);
			}
			else
				x = 0.5 * (opt->TrialPoints[iind - 1].first + opt->TrialPoints[iind].first);
			NewPoints.push_back (std::make_pair (x, 0.0));
			
			dx = pow (opt->TrialPoints[iind].first - opt->TrialPoints[iind - 1].first, 1.0 / ndim);
			if (dx < eps_stop) { // Convergence check
				stop = 1;
				break;
			}
			
#ifdef DEBUG_PRINT
			std::cout << "Inserting " << x << " to {" << opt->TrialPoints[iind - 1].first << ", " 
													  << opt->TrialPoints[iind].first << "}. Indices are: " 
													  << iind - 1 << ", " << iind << std::endl;
#endif
		}
		
		if (stop == 1)
			break;
		
// 		#pragma omp parallel for
		for (int i = 0; i < Indices.size(); ++i) {
			NewPoints[i].second = Func (opt->hc, NewPoints[i].first);
			if (NewPoints[i].second < f_min) {
				x_min = NewPoints[i].first;
				f_min = NewPoints[i].second;
			}
		}
		
#ifdef DEBUG_PRINT
		std::cout << "NewPoints:" << std::endl;
		for (auto& x: NewPoints)
			std::cout << x.first << " / " << x.second << ", ";
		std::cout << std::endl;
#endif
		
		InsertMultiple(opt->TrialPoints, NewPoints, Indices);
		nvmax = opt->TrialPoints.size();
		
#ifdef DEBUG_PRINT
		std::cout << "TrialPoints:" << std::endl;
		for (auto& x: opt->TrialPoints)
			std::cout << x.first << " / " << x.second << ", ";
		std::cout << std::endl;
#endif
	}
	STATS.close ();
	
	opt->hc->DrawHilbert (opt->TrialPoints, nvmax);
	opt->hc->Double2DoubleAxes(xN_min, x_min, 1, 1);
	
	if (iiter == niter)
		std::cout << "The method did not converge for " << iiter << " iterations" << std::endl;
	else
		std::cout << "The method converged for " << iiter << " iterations" << std::endl;
	std::cout << "Global minimum at point {" << xN_min[0];
	for (int k = 1; k < ndim; ++k)
		std::cout << ", " << xN_min[k];
	std::cout << "} is equal to " << f_min << std::endl << "Total points = " << nvmax << std::endl;
	
	time_total = clock() - time_total;
	
	std::cout << "Total time = " << (double)time_total / CLOCKS_PER_SEC << std::endl;
	delete opt;
	
	return 0;
}
