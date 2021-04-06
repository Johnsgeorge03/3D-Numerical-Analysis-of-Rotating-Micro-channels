#include "solution.hpp"

Solution::Solution():
dt(0.001), R(1.0), visc(1e-03), density(1000.0), omega(100.0),
nsteps(1), maxit(2), Imonitor(1), Jmonitor(1), Kmonitor(1),
URFUVel(0.8), URFVVel(0.8), URFWVel(0.8), URFPressure(0.2),alfa(0.92)
{
}


Solution::~Solution()
{
}
