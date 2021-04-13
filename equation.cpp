#include "equation.hpp"

Equation::Equation(const FiniteMatrix::finiteMat& Fmatrix)
:value(0.0), Residual(0.0), RSM(0.0), RESOR(4), URF(0.8), EqnName("U-momentum"), SOR(0.2),
APInitial(Fmatrix),
APNot(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
AP(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
AE(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
AW(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
AS(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
AN(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
AB(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
AT(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
rAP(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
APU(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
SP(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
SF(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
sourceInitial(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
sourceFinal(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
sourceRelaxed(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
NI(APInitial.size()),
NJ(APInitial[0].size()),
NK(APInitial[0][0].size()),
NIM(NI-1),
NJM(NJ-1),
NKM(NK-1),
Literations(5),
UE(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
UN(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
LW(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
LS(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
LPR(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
RES(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size())))
{
	forAllInternal(AP)
	{
		AE[i][j][k].value = 0.0;
		AW[i][j][k].value = 0.0;
		AN[i][j][k].value = 0.0;
		AS[i][j][k].value = 0.0;
		AT[i][j][k].value = 0.0;
		AB[i][j][k].value = 0.0;
		AP[i][j][k].value = 0.0;
		APNot[i][j][k].value = 0.0;
		SP[i][j][k].value = 0.0;
		SF[i][j][k].value = 0.0;
	}

	forAllInternal(AP)
	{
		AE[i][j][k].value	= APInitial[i][j][k].ae;
		AW[i][j][k].value 	= APInitial[i][j][k].aw;
		AN[i][j][k].value 	= APInitial[i][j][k].an;
		AS[i][j][k].value 	= APInitial[i][j][k].as;
		AT[i][j][k].value 	= APInitial[i][j][k].at;
		AB[i][j][k].value 	= APInitial[i][j][k].ab;
		APNot[i][j][k].value 	= APInitial[i][j][k].apnot;
		SP[i][j][k].value 	= APInitial[i][j][k].sp;
		SF[i][j][k].value 	= APInitial[i][j][k].sf;
		// let APInitial hold all these values initially, but eventually the values are sepated 
		// into different arrays having same link coefficient name
		
	}
	
}

void Equation::assembleEquation(Fields::vec3dField& massFluxEast, Fields::vec3dField& massFluxNorth, 
				Fields::vec3dField& massFluxTop, Fields::vec3dField& Ucell, Fields::vec3dField& Vcell,
				Fields::vec3dField& Wcell,int& iter)
{
// only for the internal cvs and not the boundary and boundary neighbours
	for(unsigned int i = 2; i<AP.size() - 2; i++)
	{
		for(unsigned int j = 2 ; j<AP[i].size() - 2; j++)
		{
			for(unsigned int k = 2; k< AP[i][j].size() - 2 ; k++)
			{
				AP[i][j][k].value = APNot[i][j][k].value + AE[i][j][k].value + AW[i][j][k].value
				    + AS[i][j][k].value + AN[i][j][k].value + AB[i][j][k].value 
				    + AT[i][j][k].value
				    + massFluxEast[i][j][k].value - massFluxEast[i-1][j][k].value
				    + massFluxNorth[i][j][k].value - massFluxNorth[i][j-1][k].value
				    + massFluxTop[i][j][k].value - massFluxTop[i][j][k-1].value;

				rAP[i][j][k].value = 1.0/AP[i][j][k].value; // reciprocal of AP
				

				// for debugging only ** not calculation****
				if(iter >= 0)
				{
					APU[i][j][k].value = (AP[i][j][k].value - (AE[i][j][k].value + AW[i][j][k].value 
								+ AS[i][j][k].value + AN[i][j][k].value + AB[i][j][k].value
								+ AT[i][j][k].value));
 				}
			}
		}
	}
	//west

	for(unsigned int i = 1; i<2;i++)
	{
		for(unsigned int j = 1; j<AP[0].size() - 1; j++)
		{
			for(unsigned int k = 1; k<AP[0][0].size() - 1; k++)
			{
			AP[i][j][k].value = APNot[i][j][k].value + AE[i][j][k].value + AW[i][j][k].value
				+ AS[i][j][k].value + AN[i][j][k].value + AB[i][j][k].value 
				+ AT[i][j][k].value
				+ massFluxEast[i][j][k].value - massFluxEast[i-1][j][k].value
				+ massFluxNorth[i][j][k].value - massFluxNorth[i][j-1][k].value
				+ massFluxTop[i][j][k].value - massFluxTop[i][j][k-1].value
				+ 8.0*Ucell[i][j][k].visc*Ucell[i][j][k].Se/(3.0*Ucell[i][j][k].DXPtoE);

			rAP[i][j][k].value = 1.0/AP[i][j][k].value;
			}
		}
	}

	// east
	for(unsigned int i = AP.size() - 2; i < AP.size() -1; i++)
	{
		for(unsigned int j = 1; j<AP[0].size() - 1; j++)
		{
			for(unsigned int k = 1; k<AP[0][0].size() - 1; k++)
			{
			AP[i][j][k].value = APNot[i][j][k].value + AE[i][j][k].value + AW[i][j][k].value
				+ AS[i][j][k].value + AN[i][j][k].value + AB[i][j][k].value 
				+ AT[i][j][k].value
				+ massFluxEast[i][j][k].value - massFluxEast[i-1][j][k].value
				+ massFluxNorth[i][j][k].value - massFluxNorth[i][j-1][k].value
				+ massFluxTop[i][j][k].value - massFluxTop[i][j][k-1].value
				+ 8.0*Ucell[i][j][k].visc*Ucell[i][j][k].Se/(3.0*Ucell[i-1][j][k].DXPtoE);

			rAP[i][j][k].value = 1.0/AP[i][j][k].value;

			}
		}
	}


	// south

	for(unsigned int i = 1; i < AP.size() - 1; i++)
	{
		for(unsigned int j = 1; j<2; j++)
		{
			for( unsigned int k = 1; k< AP[0][0].size() - 1; k++)
			{
			AP[i][j][k].value = APNot[i][j][k].value + AE[i][j][k].value + AW[i][j][k].value
				+ AS[i][j][k].value + AN[i][j][k].value + AB[i][j][k].value 
				+ AT[i][j][k].value
				+ massFluxEast[i][j][k].value - massFluxEast[i-1][j][k].value
				+ massFluxNorth[i][j][k].value - massFluxNorth[i][j-1][k].value
				+ massFluxTop[i][j][k].value - massFluxTop[i][j][k-1].value
				+ 8.0*Vcell[i][j][k].visc*Vcell[i][j][k].Sn/(3.0*Vcell[i][j][k].DYPtoN);

			rAP[i][j][k].value = 1.0/AP[i][j][k].value;
			}
		}	
	}

	//north

	for(unsigned int i = 1; i < AP.size() -1 ; i++)
	{
		for(unsigned int j = AP[0].size() - 2; j < AP[0].size() - 1; j++)
		{
			for(unsigned int k = 1; k < AP[0][0].size() - 1; k++)
			{
			AP[i][j][k].value = APNot[i][j][k].value + AE[i][j][k].value + AW[i][j][k].value
				+ AS[i][j][k].value + AN[i][j][k].value + AB[i][j][k].value 
				+ AT[i][j][k].value
				+ massFluxEast[i][j][k].value - massFluxEast[i-1][j][k].value
				+ massFluxNorth[i][j][k].value - massFluxNorth[i][j-1][k].value
				+ massFluxTop[i][j][k].value - massFluxTop[i][j][k-1].value
				+ 8.0*Vcell[i][j][k].visc*Vcell[i][j][k].Sn/(3.0*Vcell[i-1][j][k].DYPtoN);

			rAP[i][j][k].value = 1.0/AP[i][j][k].value;
			}
		}
	}

	//bottom
	for(unsigned int  i = 1; i < AP.size() - 1; i++)
	{
		for(unsigned int j = 1; j<AP[0].size() - 1; j++)
		{
			for(unsigned int k = 1 ;k<2; k++)
			{
			AP[i][j][k].value = APNot[i][j][k].value + AE[i][j][k].value + AW[i][j][k].value
				+ AS[i][j][k].value + AN[i][j][k].value + AB[i][j][k].value 
				+ AT[i][j][k].value
				+ massFluxEast[i][j][k].value - massFluxEast[i-1][j][k].value
				+ massFluxNorth[i][j][k].value - massFluxNorth[i][j-1][k].value
				+ massFluxTop[i][j][k].value
				+ 8.0*Vcell[i][j][k].visc*Vcell[i][j][k].Sn/(3.0*Vcell[i-1][j][k].DYPtoN);

			rAP[i][j][k].value = 1.0/AP[i][j][k].value;
			}
		}
	}

	//top
	for(unsigned int i = 1; i<AP.size() - 1; i++)
	{
		for(unsigned int j = 1; j<AP[0].size() - 1; j++)
		{
			for(unsigned int k = AP[0][0].size() - 2; k < AP[0][0].size() - 1; k++)
			{
			AP[i][j][k].value = APNot[i][j][k].value + AE[i][j][k].value + AW[i][j][k].value
				+ AS[i][j][k].value + AN[i][j][k].value + AB[i][j][k].value 
				+ AT[i][j][k].value
				+ massFluxEast[i][j][k].value - massFluxEast[i-1][j][k].value
				+ massFluxNorth[i][j][k].value - massFluxNorth[i][j-1][k].value
			- massFluxTop[i][j][k-1].value + Wcell[i][j][k].density*Wcell[i][j][k].St*Wcell[i][j][k].value; // no flux fromoutlet
				

			rAP[i][j][k].value = 1.0/AP[i][j][k].value;
		
			}
		}
	}
}




void Equation::assemblePressureEquation()
{
	for(unsigned int i = 1; i<AP.size() - 1; i++)
	{
		for(unsigned int j = 1; j<AP[0].size() - 1; j++)
		{
			for(unsigned int k = 1; k< AP[0][0].size() - 1; k++)
			{
				AP[i][j][k].value = AE[i][j][k].value + AW[i][j][k].value + AS[i][j][k].value
						+ AN[i][j][k].value + AT[i][j][k].value + AB[i][j][k].value;

				rAP[i][j][k].value = 1.0/AP[i][j][k].value;

			}
		}
	}

}


//Gauss-Siedel Solver

Fields::vec3dField Equation::solveVelocity(Fields::vec3dField& phi, unsigned int& iterations)
{
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), Fields::vec1dField(phi[0][0].size())));
	Fields::vec3dField phiold (phi.size(), Fields::vec2dField(phi[0].size(), Fields::vec1dField(phi[0][0].size())));
	phiold = phi;
	for(unsigned int iter = 1; iter<= iterations; iter++)
	{
		for(unsigned int i = 1; i< phi.size() - 1; i++)
		{
			for(unsigned int j = 1 ; j<phi[0].size() - 1; j++)
			{
				for(unsigned int k = 1; k<phi[0][0].size() - 1; k++)
				{
					phi[i][j][k].value = rAP[i][j][k].value*(AW[i][j][k].value*phi[i-1][j][k].value
					+ AE[i][j][k].value*phi[i+1][j][k].value + AS[i][j][k].value*phi[i][j-1][k].value
					+ AN[i][j][k].value*phi[i][j+1][k].value + AB[i][j][k].value*phi[i][j][k-1].value
					+ AT[i][j][k].value*phi[i][j][k+1].value + APNot[i][j][k].value*phiold[i][j][k].value
					+ SP[i][j][k].value + SF[i][j][k].value);
				}
			}

		}
	}
	
		
	forAll(phi)
	{
		phitemp[i][j][k].value = phi[i][j][k].value;
		
	}
	
	
	return phitemp;
}




// pass pressure
Fields::vec3dField Equation::solvePressure(Fields::vec3dField& phi, unsigned int& iterations)
{
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), Fields::vec1dField(phi[0][0].size())));

	for(unsigned int iter = 1; iter<= iterations; iter++)
	{
		for(unsigned int i = 1; i< phi.size() - 1; i++)
		{
			for(unsigned int j = 1 ; j<phi[0].size() - 1; j++)
			{
				for(unsigned int k = 1; k<phi[0][0].size() - 1; k++)
				{
					phi[i][j][k].value = rAP[i][j][k].value*(AW[i][j][k].value*phi[i-1][j][k].value
					+ AE[i][j][k].value*phi[i+1][j][k].value + AS[i][j][k].value*phi[i][j-1][k].value
					+ AN[i][j][k].value*phi[i][j+1][k].value + AB[i][j][k].value*phi[i][j][k-1].value
					+ AT[i][j][k].value*phi[i][j][k+1].value + SP[i][j][k].value);
				}
			}

		}
	}			
	
	forAll(phi)
	{
		phitemp[i][j][k].value = phi[i][j][k].value;
	}
	return phitemp;


}





Fields::vec3dField Equation::momentumInterpolation(Fields::vec3dField& vel, Fields::vec3dField& pressure, 
							int& direction)
{
		
	if(direction == 1) // U vel interpolation
	{
		Fields::vec3dField wallvel (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

		forAllInternal(wallvel)
		{
			double dP = vel[i][j][k].Se*rAP[i][j][k].value;
			double dE = vel[i+1][j][k].Se*rAP[i+1][j][k].value;
			double de = (dE + dP)/2.0;
			wallvel[i][j][k].de = de;
		//cout<< "Se " <<vel[i][j][k].Se<<" reci "<<rAP[i][j][k].value<<" dP = "<<dP<<" dE = "<<dE<<" de = "<<de<<endl;
			wallvel[i][j][k].value = (vel[i][j][k].value + vel[i+1][j][k].value)/2.0 - (SP[i][j][k].value
			+ SP[i+1][j][k].value)/2.0 + de*(pressure[i][j][k].value - pressure[i+1][j][k].value);
		}

		forWestBoundary(wallvel)
		{
			wallvel[i][j][k].value = vel[i][j][k].value;
		}	

		forEastBoundary(wallvel)
		{
			wallvel[i][j][k].value = vel[i+1][j][k].value;
		
		}
		return wallvel;
	}

	else if (direction == 2) // V vel interpolation
	{
		Fields::vec3dField wallvel (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));


		forAllInternal(wallvel)
		{
			double dP = vel[i][j][k].Sn*rAP[i][j][k].value;
			double dN = vel[i][j+1][k].Sn*rAP[i][j+1][k].value;
			double dn = (dN + dP)/2.0;
			wallvel[i][j][k].dn = dn;

			wallvel[i][j][k].value = (vel[i][j][k].value + vel[i][j+1][k].value)/2.0 - (SP[i][j][k].value
			+ SP[i][j+1][k].value)/2.0 + dn*(pressure[i][j][k].value - pressure[i][j+1][k].value);


		}

		
		forSouthBoundary(wallvel)
		{
			wallvel[i][j][k].value = vel[i][j][k].value;
		}
		
		forNorthBoundary(wallvel)
		{
		
			wallvel[i][j][k].value = vel[i][j+1][k].value;
		}

		return wallvel;


	}

	else if (direction == 3) // W vel interpolation
	{
		Fields::vec3dField wallvel (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));


		forAllInternal(wallvel)
		{
			double dP = vel[i][j][k].St*rAP[i][j][k].value;
			double dT = vel[i][j][k+1].St*rAP[i][j][k+1].value;
			double dt = (dT + dP)/2.0;
			wallvel[i][j][k].dt = dt;

			wallvel[i][j][k].value = (vel[i][j][k].value + vel[i][j][k+1].value)/2.0 - (SP[i][j][k].value
			+ SP[i][j][k+1].value)/2.0 + dt*(pressure[i][j][k].value - pressure[i][j][k+1].value);
	
		}

		forBottomBoundary(wallvel)
		{
			wallvel[i][j][k].value = vel[i][j][k].value;
		}
			


		forTopBoundary(wallvel)
		{
			wallvel[i][j][k].value = vel[i][j][k+1].value;
			wallvel[i][j][k].dt = vel[i][j][k].St*rAP[i][j][k].value/2.0;
		}

		return wallvel;
	}
	
}








//For the calculation of massFlux


Fields::vec3dField Equation::massFluxCalculation(Fields::vec3dField& wallvel, 
						Fields::vec3dField& nodevel, int& direction)
{
	
	
	if(direction == 1) // east direction mass flux
	{
		Fields::vec3dField massFlux (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
	

		double Se = nodevel[1][1][1].Se;
		double density = nodevel[1][1][1].density;
		for(unsigned int i = 0; i<massFlux.size(); i++)
		{
			for(unsigned int j = 1; j< massFlux[0].size() - 1 ; j++)
			{
				for(unsigned int k = 1; k < massFlux[0][0].size() - 1; k++)
				{
					massFlux[i][j][k].value = density*Se*wallvel[i][j][k].value;	
				}
			}
		}
		return massFlux;
	}

	else if (direction == 2) // north direction mass flux
	{
		Fields::vec3dField massFlux (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
	

		double Sn = nodevel[1][1][1].Sn;
		double density = nodevel[1][1][1].density;

		for(unsigned int i = 1; i<massFlux.size() - 1; i++)
		{
			for(unsigned int j = 0; j< massFlux[0].size(); j++)
			{
				for(unsigned int k = 1; k < massFlux[0][0].size() - 1; k++)
				{
					massFlux[i][j][k].value = density*Sn*wallvel[i][j][k].value;
	
				}
			}
		}
		return massFlux;

	}

	else if (direction == 3) // top direction mass flux
	{
		Fields::vec3dField massFlux (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));
	

		double St = nodevel[1][1][1].St;
		double density = nodevel[1][1][1].density;
		for(unsigned int i = 1; i<massFlux.size() - 1; i++)
		{
			for(unsigned int j = 1; j< massFlux[0].size() - 1; j++)
			{
				for(unsigned int k = 0; k < massFlux[0][0].size(); k++)
				{
					massFlux[i][j][k].value = density*St*wallvel[i][j][k].value;	
	
				}
			}
		}
		return massFlux;
	}

}







// For velocity correction



Fields::vec3dField Equation::faceVelocityCorrection(Fields::vec3dField& PCorr, Fields::vec3dField& UWall, 
						Fields::vec3dField& VWall, Fields::vec3dField& WWall, int& direction)
{
	if(direction == 1) // Uwall vel correction
	{
		Fields::vec3dField wallvel (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

		forAllInternal(wallvel)
		{
			wallvel[i][j][k].value = UWall[i][j][k].de*(PCorr[i][j][k].value - PCorr[i+1][j][k].value);
		}

		return wallvel;
	}

	else if (direction == 2) // Vwall vel correction
	{
		Fields::vec3dField wallvel (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));


		forAllInternal(wallvel)
		{
			wallvel[i][j][k].value = VWall[i][j][k].dn*(PCorr[i][j][k].value - PCorr[i][j+1][k].value);

		}

		return wallvel;


	}

	else if (direction == 3) // Wwall vel correction
	{
		Fields::vec3dField wallvel (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));


		forAllInternal(wallvel)
		{

			wallvel[i][j][k].value = WWall[i][j][k].dt*(PCorr[i][j][k].value - PCorr[i][j][k+1].value);
	
		}
		
		forTopBoundary(wallvel)
		{
			wallvel[i][j][k].value = WWall[i][j][k].dt*(PCorr[i][j][k].value - PCorr[i][j][k+1].value);
		}


		return wallvel;
	}

}



// For the correction of massflux




Fields::vec3dField Equation::massFluxCorrection(Fields::vec3dField& PCorr, Fields::vec3dField& UWCorr, 
						Fields::vec3dField& VWCorr, Fields::vec3dField& WWCorr, int& direction)
{
	if(direction == 1) // east direction mass flux correction
	{
		Fields::vec3dField massFlux (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
	
		forAllInternal(massFlux)
		{
			massFlux[i][j][k].value = PCorr[i][j][k].Se*PCorr[i][j][k].density*UWCorr[i][j][k].value;
		}
		
		return massFlux;
	}

	else if (direction == 2) // north direction mass flux
	{
		Fields::vec3dField massFlux (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
	
		forAllInternal(massFlux)
		{
			massFlux[i][j][k].value = PCorr[i][j][k].Sn*PCorr[i][j][k].density*VWCorr[i][j][k].value;
		}
		
		return massFlux;

	}

	else if (direction == 3) // top direction mass flux
	{
		Fields::vec3dField massFlux (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));
		
		forAllInternal(massFlux)
		{
			massFlux[i][j][k].value = PCorr[i][j][k].St*PCorr[i][j][k].density*WWCorr[i][j][k].value;
		}
		forTopBoundary(massFlux)
		{
			massFlux[i][j][k].value = PCorr[i][j][k].St*PCorr[i][j][k].density*WWCorr[i][j][k].value;
		}

		return massFlux;

	}


}


Fields::vec3dField Equation::cellVelocityCorrection(Fields::vec3dField& PCorr, Fields::vec3dField& Ucentre, 
						Fields::vec3dField& Vcentre, Fields::vec3dField& Wcentre, int& direction)
{
	if(direction == 1) // U vel correction
	{
		Fields::vec3dField Ucell (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

		forAllInternal(Ucell)
		{
			Ucell[i][j][k].value = rAP[i][j][k].value*Ucentre[i][j][k].Se
						*(PCorr[i-1][j][k].value - PCorr[i+1][j][k].value);
		}

		return Ucell;
	}

	else if (direction == 2) // V vel correction
	{
		Fields::vec3dField Vcell (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));


		forAllInternal(Vcell)
		{
			Vcell[i][j][k].value = rAP[i][j][k].value*Vcentre[i][j][k].Sn
						*(PCorr[i][j-1][k].value - PCorr[i][j+1][k].value);

		}

		return Vcell;


	}

	else if (direction == 3) // W vel correction
	{
		Fields::vec3dField Wcell (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));


		forAllInternal(Wcell)
		{

			Wcell[i][j][k].value = rAP[i][j][k].value*Wcentre[i][j][k].St
						*(PCorr[i][j][k-1].value - PCorr[i][j][k+1].value);
	
		}

		return Wcell;
	}

}





Equation::~Equation()
{


}
