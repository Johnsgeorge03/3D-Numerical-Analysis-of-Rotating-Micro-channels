#include "mesh.hpp"

//NX      -  no of cells in x-direction
//lengthX -  length of the micro-channel in x-direction
//dx      -  spacing between nodes in the x direction
//FX      -  interpolation factor in the x.
//X       -  the cell face value (coordinates)
//XC      -  cell centre values (coordinates)
//beginX  -  the starting coordinate of the microchannel
//NI      -  the total no of nodes including boundary
//initializing constructor
Mesh::Mesh(double& beginX_, double& beginY_, double& beginZ_, int& NX_, int& NY_, int& NZ_, 
		double& lengthX_, double& lengthY_, double& lengthZ_):

beginX(beginX_), beginY(beginY_), beginZ(beginZ_),
NX(NX_), NY(NY_), NZ(NZ_), 
lengthX(lengthX_), lengthY(lengthY_), lengthZ(lengthZ_),
X(NX_ + 2, 0.0), Y(NY_ + 2, 0.0), Z(NZ_ + 2, 0.0), 
XC(NX_ + 2, 0.0), YC(NY_ + 2, 0.0), ZC(NZ_ + 2, 0.0),
FX(NX_ + 2, 0.0), FY(NY_ + 2, 0.0), FZ(NZ_ + 2, 0.0)
{
    NI  = NX + 2;
    NJ  = NY + 2;
    NK  = NZ + 2;
    NIM = NI - 1;
    NJM = NJ - 1;
    NKM = NK - 1;

    dx = lengthX / NX;
    dy = lengthY / NY;
    dz = lengthZ / NZ;

    setX(X);
    setY(Y);
    setZ(Z);
    setXC(X, XC);
    setYC(Y, YC);
    setZC(Z, ZC);
    setFX(X, XC, FX);
    setFY(Y, YC, FY);
    setFZ(Z, ZC, FZ);
}


Mesh::~Mesh()
{

}


void Mesh::setX(vector<double>& vecX)
{

    vecX[0] = beginX;
    for(unsigned int i = 1 ; i<vecX.size(); i++)
    {
        vecX[i] =  vecX[i-1] + dx;
    }
    vecX[vecX.size() - 1] =  vecX[vecX.size() - 2];

}


void Mesh::setY(vector<double>& vecY)
{

    vecY[0] = beginY;
    for(unsigned int i = 1 ; i<vecY.size(); i++)
    {
        vecY[i] =  vecY[i-1] + dy;
    }
    vecY[vecY.size() - 1] =  vecY[vecY.size() - 2];

}



void Mesh::setZ(vector<double>& vecZ)
{

    vecZ[0] = beginZ;
    for(unsigned int i = 1 ; i<vecZ.size(); i++)
    {
        vecZ[i] =  vecZ[i-1] + dz;
    }
    vecZ[vecZ.size() - 1] =  vecZ[vecZ.size() - 2];

}




void Mesh::setXC(vector<double>& vecX, vector<double>& vecXC)
{

    vecXC[0] = beginX;
    for(unsigned int i =  1; i < vecXC.size(); i++)
    {
        vecXC[i] = (vecX[i] + vecX[i-1])*0.5;
    }
    vecXC[vecXC.size() - 1] =  vecX[vecX.size() - 2];

}




void Mesh::setYC(vector<double>& vecY, vector<double>& vecYC)
{

    vecYC[0] = beginY;
    for(unsigned int i =  1; i < vecYC.size(); i++)
    {
        vecYC[i] = (vecY[i] + vecY[i-1])*0.5;
    }
    vecYC[vecYC.size() - 1] =  vecY[vecY.size() - 2];

}




void Mesh::setZC(vector<double>& vecZ, vector<double>& vecZC)
{

    vecZC[0] = beginZ;
    for(unsigned int i =  1; i < vecZC.size(); i++)
    {
        vecZC[i] = (vecZ[i] + vecZ[i-1])*0.5;
    }
    vecZC[vecZC.size() - 1] =  vecZ[vecZ.size() - 2];

}





void Mesh::setFX(vector<double>& vecX, vector<double>& vecXC, vector<double>& vecFX)
{
    for(unsigned int i = 0; i< vecX.size(); i++)
    {
        vecFX[i] = (vecX[i] - vecXC[i])/(vecXC[i+1] - vecXC[i]);
    }
}




void Mesh::setFY(vector<double>& vecY, vector<double>& vecYC, vector<double>& vecFY)
{
    for(unsigned int i = 0; i< vecY.size(); i++)
    {
        vecFY[i] = (vecY[i] - vecYC[i])/(vecYC[i+1] - vecYC[i]);
    }

}





void Mesh::setFZ(vector<double>& vecZ, vector<double>& vecZC, vector<double>& vecFZ)
{
    for(unsigned int i = 0; i< vecZ.size(); i++)
    {
        vecFZ[i] = (vecZ[i] - vecZC[i])/(vecZC[i+1] - vecZC[i]);
    }

}




int Mesh::getNI()
{
    return NI;
}




int Mesh::getNJ()
{
    return NJ;
}




int Mesh::getNK()
{
    return NK;
}


