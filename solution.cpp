#include "solution.hpp"

Solution::Solution():
dt(0.001), R(1.0), visc(0.001), density(1000.0), omega(0.0),
nsteps(1), maxit(50000), Imonitor(1), Jmonitor(1), Kmonitor(1),
URFUVel(0.8), URFVVel(0.8), URFWVel(0.8), URFPressure(0.6),alfa(0.9) // alfa for sip-solver
{
}


Solution::~Solution()
{
}
