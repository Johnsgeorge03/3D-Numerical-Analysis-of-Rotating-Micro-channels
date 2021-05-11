#ifndef EQUATION_H
#define EQUATION_H

#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <stdlib.h>
#include "mesh.hpp"
#include "fields.hpp"
#include "finitematrix.hpp"

class Equation
{
public:
	Equation(const FiniteMatrix::finiteMat&);
	virtual ~Equation();

	typedef vector<FiniteMatrix> Svector1d;
	typedef FiniteMatrix::finiteMat Svector; //same as that of finiteMat, just a name change
	void assembleEquation(Fields::vec3dField&, Fields::vec3dField&, Fields::vec3dField&,
			Fields::vec3dField&, Fields::vec3dField&, Fields::vec3dField&,int&);

	void relax(Fields::vec3dField&, Fields::vec3dField&);

	void assemblePressureCorrectionEquation(FiniteMatrix::finiteMat&, Fields::vec3dField&);

	//solvers
	Fields::vec3dField solveVelocity(Fields::vec3dField&, unsigned int& );
	Fields::vec3dField solvePressure(Fields::vec3dField&, unsigned int& );
	Fields::vec3dField sipSolver(Fields::vec3dField&, Fields::vec3dField&,
			            Solution&, unsigned int&, double&);
	Fields::vec3dField sipSolverPressure(Fields::vec3dField&, Solution&, unsigned int&,
						double&);
	Fields::vec3dField CGSTAB(Fields::vec3dField&, Solution&, unsigned int&, double&);

	Fields::vec3dField CGSTABvel(Fields::vec3dField&, Fields::vec3dField&, Solution&, 
						unsigned int&, double&);

	Fields::vec3dField ICCG(Fields::vec3dField&, Solution&, unsigned int&, double&);

	
	Fields::vec3dField momentumInterpolation(Fields::vec3dField&, Fields::vec3dField&, 
		Fields::vec3dField&, Fields::vec3dField&, Fields::vec3dField&, Fields::vec3dField&,
		Fields::vec3dField&, Fields::vec3dField&, Solution&, int&);

	Fields::vec3dField massFluxCalculation(Fields::vec3dField&, Fields::vec3dField&, int&);

	Fields::vec3dField faceVelocityCorrection(Fields::vec3dField&, Fields::vec3dField&, 
			Fields::vec3dField&, Fields::vec3dField&, int&);
	Fields::vec3dField massFluxCorrection(Fields::vec3dField&, Fields::vec3dField&,
			Fields::vec3dField&, Fields::vec3dField&, int&);
	Fields::vec3dField cellVelocityCorrection(Fields::vec3dField&, int&);
	double value;
	double Residual;
	double R2;
	double R2norm;
	double RESOR;

	double URF;
	std::string EqnName;
	double SOR;

	FiniteMatrix::finiteMat APInitial;
	FiniteMatrix::finiteMat APNot;
        FiniteMatrix::finiteMat AP;
        FiniteMatrix::finiteMat AE;
        FiniteMatrix::finiteMat AW;
        FiniteMatrix::finiteMat AS;
        FiniteMatrix::finiteMat AN;
        FiniteMatrix::finiteMat AB;
        FiniteMatrix::finiteMat AT;
        FiniteMatrix::finiteMat rAP;  //reciprocal of AP
        FiniteMatrix::finiteMat APU;
	FiniteMatrix::finiteMat SP;
	FiniteMatrix::finiteMat SF;
        FiniteMatrix::finiteMat sourceInitial;
        FiniteMatrix::finiteMat sourceFinal;
	FiniteMatrix::finiteMat sourceRelaxed;


	int NI, NJ, NK, NIM, NJM, NKM, Literations;
	Svector UE, UN, UT, LB, LW, LS, LPR, RES, AUX, DELTA;

	//for CGSTAB
	double ALF, BETO, GAM;
	Svector D, RESO, PK, UK, ZK, VK;

	// for testing
	FiniteMatrix testmat;
	Fields testfield;
};

#endif
