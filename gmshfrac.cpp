#include <gmsh.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "gmshfrac.h"

int gmshFracGrid (double * xy, const char * fname)
{
	double beg[2] = {xy[0], xy[1]}, end[2] = {xy[2], xy[3]};
	double length = sqrt ((beg[0] - end[0]) * (beg[0] - end[0]) + (beg[1] - end[1]) * (beg[1] - end[1]));

	std::cout << "unit fracture is: (" << beg[0] << ", " << beg[1] << ") - (" << end[0] << ", " << end[1] << ")" << std::endl;
	std::cout << "fracture length = " << length << std::endl;

	// map from [0,1]
	double eps = 0.001;
	double lx = 0.15 + eps;
	double rx = 1.85 - eps;
	double dx = rx - lx;
	double ly = 0.0 + eps;
	double ry = 1.0 - eps;
	double dy = ry - ly;

	int npw = 4;
	int nplrs = 21;
	int npbts = 38;
	int npf = 41;

	double begs[2] = {beg[0] * dx + lx, beg[1] * dy + ly};
	double ends[2] = {end[0] * dx + lx, end[1] * dy + ly};

	gmsh::initialize();
	gmsh::model::add("frac");

	gmsh::model::geo::addPoint(0.15, 0, 0, 1, 1);
	gmsh::model::geo::addPoint(0.15, 1, 0, 1, 2);
	gmsh::model::geo::addPoint(1.85, 1, 0, 1, 3);
	gmsh::model::geo::addPoint(1.85, 0, 0, 1, 4);

	gmsh::model::geo::addPoint(begs[0], begs[1], 0.0, 1, 5); 
	gmsh::model::geo::addPoint(ends[0], ends[1], 0.0, 1, 6);

	gmsh::model::geo::addPoint(0, 0, 0, 1, 9);
	gmsh::model::geo::addPoint(0, 1, 0, 1, 10);

	gmsh::model::geo::addPoint(2, 0, 0, 1, 11);
	gmsh::model::geo::addPoint(2, 1, 0, 1, 12);

	gmsh::model::geo::addLine(1, 2, 1);
	gmsh::model::geo::addLine(2, 3, 2);
	gmsh::model::geo::addLine(3, 4, 3);
	gmsh::model::geo::addLine(4, 1, 4);

	gmsh::model::geo::addLine(5, 6, 5);

	gmsh::model::geo::addLine(2, 10, 7);
	gmsh::model::geo::addLine(10, 9, 8);
	gmsh::model::geo::addLine(9, 1, 9);

	gmsh::model::geo::addLine(4, 11, 10);
	gmsh::model::geo::addLine(11, 12, 11);
	gmsh::model::geo::addLine(12, 3, 12);

	gmsh::model::geo::addCurveLoop({2, 3, 4, 1}, 13);
	gmsh::model::geo::addCurveLoop({7, 8, 9, 1}, 14);
	gmsh::model::geo::addCurveLoop({10, 11, 12, 3}, 15);

	gmsh::model::geo::addPlaneSurface({13}, 1);
	gmsh::model::geo::addPlaneSurface({14}, 2);
	gmsh::model::geo::addPlaneSurface({15}, 3);

	gmsh::model::geo::synchronize();
	gmsh::model::mesh::embed(1, {5}, 2, 1);

	gmsh::model::geo::mesh::setTransfiniteCurve(1, nplrs);
	gmsh::model::geo::mesh::setTransfiniteCurve(3, nplrs);
	gmsh::model::geo::mesh::setTransfiniteCurve(8, nplrs);
	gmsh::model::geo::mesh::setTransfiniteCurve(11, nplrs);

	gmsh::model::geo::mesh::setTransfiniteCurve(2, npbts);
	gmsh::model::geo::mesh::setTransfiniteCurve(4, npbts);

	gmsh::model::geo::mesh::setTransfiniteCurve(7, npw);
	gmsh::model::geo::mesh::setTransfiniteCurve(9, npw);
	gmsh::model::geo::mesh::setTransfiniteCurve(10, npw);
	gmsh::model::geo::mesh::setTransfiniteCurve(12, npw);

	gmsh::model::geo::mesh::setTransfiniteCurve(7, 4);

	gmsh::model::geo::mesh::setTransfiniteSurface(2, {10, 2, 1, 9});
	gmsh::model::geo::mesh::setTransfiniteSurface(3, {3, 12, 11, 4});

	gmsh::model::geo::mesh::setRecombine(2, 2);
	gmsh::model::geo::mesh::setRecombine(2, 3);

	gmsh::vectorpair outDimTags;
	gmsh::model::geo::extrude({{1, 5}}, 0, 0, 0.1, outDimTags, {1}, {1}, true);

	gmsh::model::geo::mesh::setTransfiniteCurve(5, npf);
	gmsh::model::geo::mesh::setTransfiniteCurve(16, npf);

	gmsh::model::geo::extrude({{2,1},{2,2},{2,3}}, 0, 0, 0.1, outDimTags, {1}, {1}, true);

	gmsh::model::geo::synchronize();
	gmsh::model::mesh::embed(1, {16}, 2, 41);

	gmsh::model::geo::synchronize();
	gmsh::model::mesh::embed(2, {19}, 3, 1);

	gmsh::model::addPhysicalGroup(2, {19}, 1000, "fracture");

	gmsh::model::geo::synchronize();

	std::vector<std::string> log;
	gmsh::logger::start();
	gmsh::model::mesh::generate(3);
	gmsh::logger::get(log);
	gmsh::logger::stop();
	std::string log_line = log.back();
	std::string del1 = "nodes ";
	std::string del2 = " elements";
	std::string str_ncells;
	str_ncells = log_line.substr (log_line.find(del1) + del1.size(), log_line.find(del2) - (log_line.find(del1) + del1.size()));
	int ncells = stoi(str_ncells);
	
	gmsh::option::setNumber("Mesh.SaveAll", 1);
	gmsh::write(fname);
	
	gmsh::finalize();
	
	return ncells;
}

void addFracMaterial (int ncells, const char * fname)
{
	std::ofstream fs;
	fs.open (fname, std::ios_base::app);

	fs << std::endl;
	fs << "SCALARS PERM double 1" << std::endl;
	fs << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < ncells; ++i)
		fs << "1.0" << std::endl;
	fs.close ();
}
