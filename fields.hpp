#ifndef FIELDS_H
#define FIELDS_H
#include "foralloperations.hpp"
#include "mesh.hpp"
#include "solution.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::string;

class Fields
{
public:
    // default constructor

    Fields();

    // overloading constructor

    Fields ( int&, int&, int&);
    virtual ~Fields();

    typedef vector<Fields> vec1dField;
    typedef vector<vec1dField> vec2dField;
    typedef vector<vec2dField> vec3dField;
   
    double value;
    int NI, NJ, NK, NIM, NJM, NKM;
    double de, dn, dt; // for massflux and pressure correction
    double X, XC, Y, YC, Z, ZC, FXE, FXP, FYN, FYP, FZT, FZP, DXPtoE, DYPtoN, DZPtoT; // from the mesh
    double Se, Sn, St, visc, density, volume; // Se- face area east
    
    void copyAllField(Fields::vec3dField&, Fields::vec3dField&);
    void getGridInfoPassed(Fields::vec3dField&, Mesh&, Solution&);
    void setVectorFieldGridFeatures();

    void copyInternalField(Fields::vec3dField&, Fields::vec3dField&);
    void initializeField(Fields::vec3dField&, double);
    void initializeInternalField(Fields::vec3dField&, double);
    void print3dmat(const Fields::vec3dField&);
    //boundary conditions which modifieds the values of fields
    void boundaryCondition(Fields::vec3dField&, string&, double);
    void linearextrapolateCondition(Fields::vec3dField&, vector<double>&, vector<double>&, 
    					vector<double>&, string&);
    void copyOutletVelocity(Fields::vec3dField&);

    //Fields::vec3dField& operator=(const Fields::vec3dField&); 
	
};

#endif
    
