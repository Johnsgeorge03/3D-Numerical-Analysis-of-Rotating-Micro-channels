#ifndef EQUATION_H
#define EQUATION_H

#include <string>
#include <vector>
#include "mesh.hpp"
#include "fields.hpp"
#include "finitematrix.hpp"

class Equation
{
public:
	Equation(const FiniteMatrix::finiteMat&);
	virtual ~Equation();

	typedef vector<FiniteMatrix> Svector1d;
	typedef FiniteMatrix::finiteMat Svector; // same as that of finiteMat, just a name change

	//void relax(Fields::vec3dField&);
	//void resetEqn();
	//void noWallShearXBoundaryCondition(Fields::vec3dField&);
	//void noWallShearYBoundaryCondition(Fields::vec3dField&);
	//void noWallShearZBoundaryCondition(Fields::vec3dField&);
	void assembleEquation(Fields::vec3dField&, Fields::vec3dField&, Fields::vec3dField&
				,int&);
	void assemblePressureEquation();

	Fields::vec3dField solveVelocity(Fields::vec3dField&, unsigned int& );
	Fields::vec3dField solvePressure(Fields::vec3dField&, unsigned int& );
	
	Fields::vec3dField momentumInterpolation(Fields::vec3dField&, Fields::vec3dField&, int&);
	Fields::vec3dField massFluxCalculation(Fields::vec3dField&, Fields::vec3dField&, int&);

	Fields::vec3dField velocityCorrection(Fields::vec3dField&, Fields::vec3dField&, 
			Fields::vec3dField&, Fields::vec3dField&, int&);
	Fields::vec3dField massFluxCorrection(Fields::vec3dField&, Fields::vec3dField&,
			Fields::vec3dField&, Fields::vec3dField&, int&);
	
	double value;
	double Residual;
	double RSM; // root-square name
	double RESOR;

	double URF, rURF;
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
        FiniteMatrix::finiteMat rAP;
        FiniteMatrix::finiteMat APU;
	FiniteMatrix::finiteMat SP;
	FiniteMatrix::finiteMat SF;
        FiniteMatrix::finiteMat sourceInitial;
        FiniteMatrix::finiteMat sourceFinal;
	FiniteMatrix::finiteMat sourceRelaxed;

private:
	int NI, NJ, NK, NIM, NJM, NKM, Literations;
	Svector UE, UN, LW, LS, LPR, RES;
};

#endif
