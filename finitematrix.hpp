#ifndef FINITEMATRIX_H
#define FINITEMATRIX_H

#include <vector>
#include <iostream>
#include "foralloperations.hpp"
#include "fields.hpp"

using namespace std;

class FiniteMatrix
{
public:
	FiniteMatrix();
	virtual ~FiniteMatrix();

	typedef vector<vector<vector<FiniteMatrix> > > finiteMat;//a 3d array of class finitematrix
	
	double value, apnot, ae, aw, an, as, at, ab, sp, sf; // the link coefficients and source term
		
	void print3dmat(finiteMat&);
	void print3dsource(finiteMat&);
	

	//operator overloading
	
	friend FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat&, 
						const FiniteMatrix::finiteMat&);
	friend FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat&, 
						const FiniteMatrix::finiteMat&);
	friend FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat&, 
						const FiniteMatrix::finiteMat&);
	friend FiniteMatrix::finiteMat operator*(const double, 
						const FiniteMatrix::finiteMat&);

};

#endif

// value = ap
// sp = pressure source
// sf = force source
// apnot = coeff of phiPnot in current timestep
