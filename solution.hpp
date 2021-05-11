#ifndef SOLUTION_H
#define SOLUTION_H

class Solution
{   
public:
    Solution();
    virtual ~Solution();
    double dt;
    double R;
    double visc;
    double density;
    double omega;
    int nsteps, maxit;
    int Imonitor, Jmonitor, Kmonitor;
    double URFUVel, URFVVel, URFWVel, URFPressure; // under-relaxation factor for U,V,W,P
    double URFU, URFV, URFW, URFP; // reciprocal of the above solvers
    double alfa; // for sip solver
};

#endif

// dt - time step size
// nsteps - no of time steps
// maxit - maximum iterations
// URFUVel - under-relaxation factor for u  velocity
// URFU - reciprocal of URFUVel
// omega - rps
