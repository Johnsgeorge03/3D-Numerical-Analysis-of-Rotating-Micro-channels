#ifndef MESH_H
#define MESH_H

#include<vector>

using std::vector;

class Mesh
{
private:
    double beginX, beginY, beginZ;
    int NX, NY, NZ, NI, NJ, NK, NIM , NJM, NKM;
    double lengthX, lengthY, lengthZ;
    double dx, dy, dz;

public:
    Mesh(double&, double&, double&, int&, int&, int&, double&, double&, double&);
    virtual ~Mesh();

    vector<double> X;
    vector<double> Y;
    vector<double> Z;
    vector<double> XC;
    vector<double> YC;
    vector<double> ZC;
    vector<double> FX;
    vector<double> FY;
    vector<double> FZ;

    void setX(vector<double>&);
    void setY(vector<double>&);
    void setZ(vector<double>&);
    void setXC(vector<double>&, vector<double>&);
    void setYC(vector<double>&, vector<double>&);
    void setZC(vector<double>&, vector<double>&);
    void setFX(vector<double>&, vector<double>&, vector<double>&);
    void setFY(vector<double>&, vector<double>&, vector<double>&);
    void setFZ(vector<double>&, vector<double>&, vector<double>&);
    int getNI();
    int getNJ();
    int getNK();

};

#endif
