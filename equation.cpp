#include "equation.hpp"

Equation::Equation(const FiniteMatrix::finiteMat& Fmatrix)
:value(0.0), Residual(0.0), R2(0.0), R2norm(0.0), RESOR(0.0), URF(0.8), EqnName("U-Eqn"), SOR(0.05),
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
UT(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
LB(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
LW(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
LS(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
LPR(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
RES(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
AUX(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
DELTA(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),

ALF(1.0),
BETO(1.0),
GAM(1.0),

D(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
RESO(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
PK(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
UK(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
ZK(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
VK(Fmatrix.size(),vector<vector<FiniteMatrix>>(Fmatrix[0].size(), vector<FiniteMatrix>(Fmatrix[0][0].size()))),
testfield(NI, NJ, NK)
{
	forAll(AP)
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

//--------------------------------- ASSEMBLE MOMENTUM EQUATION ----------------------------------------------//

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
				
				assert(std::isinf(AP[i][j][k].value)!=1);

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
				
			assert(std::isinf(AP[i][j][k].value)!=1);

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
				assert(std::isinf(AP[i][j][k].value)!=1);

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
				+ 8.0*Vcell[i][j][k].visc*Vcell[i][j][k].Sn/(3.0*Vcell[i][j][k].DYPtoN);

			if(std::isinf(AP[i][j][k].value) == 1)
			{
				cout<<" inf value in "<<i<<" "<<j<<" "<<k<<endl;
			}
			
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
				+ 8.0*Wcell[i][j][k].visc*Wcell[i][j][k].St/(3.0*Wcell[i][j][k].DZPtoT);

			if(std::isinf(AP[i][j][k].value) == 1)
			{
				cout<<" inf value in "<<i<<" "<<j<<" "<<k<<endl;
			}
	

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
			- massFluxTop[i][j][k-1].value 
			+ Wcell[i][j][k].density*Wcell[i][j][k].St*Wcell[i][j][k].value; 
			// no flux fromoutlet
				
			assert(std::isinf(AP[i][j][k].value)!= 1);
			rAP[i][j][k].value = 1.0/AP[i][j][k].value;
		
			}
		}
	}
} 

//------------------------------------ END OF MOMENTUM ASSEMBLY ----------------------------------------------//




//--------------------------------------- UNDER-RELAXATION ---------------------------------------------------//

void Equation::relax(Fields::vec3dField& vec, Fields::vec3dField& vecold)
{
	forAllInternal(AP)
	{
		AP[i][j][k].value = (1.0/URF)*AP[i][j][k].value;
		sourceRelaxed[i][j][k].value = SP[i][j][k].value + SF[i][j][k].value
						+ APNot[i][j][k].value*vecold[i][j][k].value
						+ (1.0-URF)*(AP[i][j][k].value*vec[i][j][k].value);
						
		rAP[i][j][k].value = 1/AP[i][j][k].value;
	}


}

//---------------------------------- END OF UNDER-RELAXATION -------------------------------------------------//




//----------------------------ASSEMBLE PRESSURECORRECTION EQUATION -------------------------------------------//

void Equation::assemblePressureCorrectionEquation(FiniteMatrix::finiteMat& APW, Fields::vec3dField& vec)
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

		//for cell near outlet boundary (top)

	for(unsigned int i = 1; i<AP.size()-1; i++)
	{
		for(unsigned int j = 1; j<AP[0].size() - 1; j++)
		{
			for(unsigned int  k = AP[0][0].size() - 2; k < AP[0][0].size() - 1; k++)
			{
				AP[i][j][k].value = AP[i][j][k].value 
			 + vec[i][j][k].density*vec[i][j][k].St*vec[i][j][k].St/APW[i][j][k].value;
			 
			 	rAP[i][j][k].value = 1.0/AP[i][j][k].value;
			}
		}
	
	}


}

// ---------------------------------- END OF PRESSURE CORRECTION ASSEMBLY -----------------------------------//





// -------------------------------------GAUSS-SIEDEL SOLVER FOR VELOCITY ------------------------------------//

Fields::vec3dField Equation::solveVelocity(Fields::vec3dField& phi, unsigned int& iterations)
{
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), 
					Fields::vec1dField(phi[0][0].size())));
	Fields::vec3dField phiold (phi.size(), Fields::vec2dField(phi[0].size(), 
					Fields::vec1dField(phi[0][0].size())));
	phiold = phi;
	for(unsigned int iter = 1; iter<= iterations; iter++)
	{
		for(unsigned int i = 1; i< phi.size() - 1; i++)
		{
			for(unsigned int j = 1 ; j<phi[0].size() - 1; j++)
			{
				for(unsigned int k = 1; k<phi[0][0].size() - 1; k++)
				{
					phi[i][j][k].value = rAP[i][j][k].value
					*(AW[i][j][k].value*phi[i-1][j][k].value
					+ AE[i][j][k].value*phi[i+1][j][k].value 
					+ AS[i][j][k].value*phi[i][j-1][k].value
					+ AN[i][j][k].value*phi[i][j+1][k].value 
					+ AB[i][j][k].value*phi[i][j][k-1].value
					+ AT[i][j][k].value*phi[i][j][k+1].value 
					+ APNot[i][j][k].value*phiold[i][j][k].value
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
//-----------------------------END OF GAUSS-SIEDEL SOLVER FOR VELOCITY---------------------------------------//



//------------------------------- GAUSS-SIEDEL SOLVER FOR PRESSURE ------------------------------------------//
// pass pressure
Fields::vec3dField Equation::solvePressure(Fields::vec3dField& phi, unsigned int& iterations)
{
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), Fields::vec1dField(phi[0][0].size())));


	
	for(unsigned int iter = 1; iter<= iterations; iter++)
	{
	Residual = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				RES[i][j][k].value = SP[i][j][k].value
				- (AP[i][j][k].value*phi[i][j][k].value
				- AE[i][j][k].value*phi[i+1][j][k].value 
				- AW[i][j][k].value*phi[i-1][j][k].value
				- AS[i][j][k].value*phi[i][j-1][k].value 
				- AN[i][j][k].value*phi[i][j+1][k].value
				- AT[i][j][k].value*phi[i][j][k+1].value 
				- AB[i][j][k].value*phi[i][j][k-1].value);

				Residual += std::pow(RES[i][j][k].value, 2.0);

			}

		}

	}

	R2 = std::sqrt(Residual);
	double small = 1e-20;
	if(iter == 1)
	{
		RESOR = R2;
	}
	
	R2norm = R2/(RESOR + small);
	
	cout<<EqnName<< " Inner iteration of gauss siedel Solver "<< iter+1
			<<" and R2 --> " <<R2<<" R2norm "<< R2norm <<endl;

		for(unsigned int i = 1; i< phi.size() - 1; i++)
		{
			for(unsigned int j = 1 ; j<phi[0].size() - 1; j++)
			{
				for(unsigned int k = 1; k<phi[0][0].size() - 1; k++)
				{
					phi[i][j][k].value = rAP[i][j][k].value
					*(AW[i][j][k].value*phi[i-1][j][k].value
					+ AE[i][j][k].value*phi[i+1][j][k].value 
					+ AS[i][j][k].value*phi[i][j-1][k].value
					+ AN[i][j][k].value*phi[i][j+1][k].value 
					+ AB[i][j][k].value*phi[i][j][k-1].value
					+ AT[i][j][k].value*phi[i][j][k+1].value 
					+ SP[i][j][k].value);
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

//---------------------------- END OF GAUSS - SIEDEL SOLVER FOR PRESSURE -------------------------------------//




//------------------------------------ SIP SOLVER FOR VELOCITY -----------------------------------------------//

Fields::vec3dField Equation::sipSolver(Fields::vec3dField& phi, Fields::vec3dField& phiold, Solution& sol, 
					unsigned int& iterations, double& tolerance)

{
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), 
					Fields::vec1dField(phi[0][0].size())));

	// coefficients of upper and lower triangular matrix
	
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				LB[i][j][k].value = -AB[i][j][k].value
							/( 1.0 + sol.alfa*(UN[i][j][k-1].value 
							+ UE[i][j][k-1].value));
				LW[i][j][k].value = -AW[i][j][k].value
							/( 1.0 + sol.alfa*(UN[i-1][j][k].value 
							+ UT[i-1][j][k].value));
				LS[i][j][k].value = -AS[i][j][k].value
							/( 1.0 + sol.alfa*(UE[i][j-1][k].value 
							+ UT[i][j-1][k].value));

				double H1 = sol.alfa*(LB[i][j][k].value*UN[i][j][k-1].value
							+ LW[i][j][k].value*UN[i-1][j][k].value);
				double H2 = sol.alfa*(LB[i][j][k].value*UE[i][j][k-1].value
							+ LS[i][j][k].value*UE[i][j-1][k].value);
				double H3 = sol.alfa*(LW[i][j][k].value*UT[i-1][j][k].value
							+ LS[i][j][k].value*UT[i][j-1][k].value);

 				LPR[i][j][k].value = 1.0/(AP[i][j][k].value + H1 + H2 + H3 
				- LB[i][j][k].value*UT[i][j][k-1].value- LW[i][j][k].value*UE[i-1][j][k].value 
				- LS[i][j][k].value*UN[i][j-1][k].value + 1e-30);
					// inertial damping
			

				UN[i][j][k].value = (-AN[i][j][k].value - H1)*LPR[i][j][k].value;
				UE[i][j][k].value = (-AE[i][j][k].value - H2)*LPR[i][j][k].value;
				UT[i][j][k].value = (-AT[i][j][k].value - H3)*LPR[i][j][k].value;
			}

		}

	}
	

	//calculating and over-riding the residual
	for(unsigned int L = 0; L<iterations; L++)
	{
	Residual = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				RES[i][j][k].value = sourceRelaxed[i][j][k].value
				+ APNot[i][j][k].value*phiold[i][j][k].value 
				- (AP[i][j][k].value*phi[i][j][k].value
				- AE[i][j][k].value*phi[i+1][j][k].value 
				- AW[i][j][k].value*phi[i-1][j][k].value
				- AS[i][j][k].value*phi[i][j-1][k].value 
				- AN[i][j][k].value*phi[i][j+1][k].value
				- AT[i][j][k].value*phi[i][j][k+1].value 
				- AB[i][j][k].value*phi[i][j][k-1].value);

				Residual += std::pow(RES[i][j][k].value, 2.0);

				AUX[i][j][k].value = (RES[i][j][k].value 
						- LB[i][j][k].value*AUX[i][j][k-1].value
						- LW[i][j][k].value*AUX[i-1][j][k].value 
						- LS[i][j][k].value*AUX[i][j-1][k].value)*LPR[i][j][k].value; 
						// U*delta
				
				//LUdelta = RES
				//U*delta = Linv* RES
			}

		}

	}


/*	// TESTING ONLY

	cout<<" ITERATION NO: "<<L + 1<<endl;
	cout<<EqnName<<"  - RES " <<endl;
	testmat.print3dmat(RES);
	cout<<" END RES "<<endl;
*/	// TESTING END

	R2 = std::sqrt(Residual);
	double small = 1e-30;
	if(L==0)
	{
		RESOR = R2;
	}
	
	R2norm = R2/(RESOR + small);
	
	//Back substitution and correction
	for(unsigned int k = phi[0][0].size() - 2; k >= 1; --k)
	{
		for(unsigned int i = phi.size() - 2; i >= 1; --i)
		{
			for(unsigned int j = phi[0].size() - 2; j >= 1; --j)
			{
				DELTA[i][j][k].value = AUX[i][j][k].value 
						- UN[i][j][k].value*DELTA[i][j+1][k].value
						- UE[i][j][k].value*DELTA[i+1][j][k].value
						- UT[i][j][k].value*DELTA[i][j][k+1].value; 

				phi[i][j][k].value = phi[i][j][k].value + DELTA[i][j][k].value;
			}


		}

	}
	// check convergence
	cout<<EqnName<< " Inner iteration of SIP Solver "<< L+1
				<< " and R2 --> " <<R2<< " R2norm "<< R2norm <<endl;
	if(R2norm <= SOR)
	{
		cout<<EqnName<<": Inner iteration of SIP Solver converged at iteration "<< L+1<<endl;
		break;
	}
/*
	if(R2norm <= tolerance)
	{
		cout<<EqnName<<": Inner iteration of SIP Solver converged at iteration "<< L+1<<endl;
		break;
	}

	else if(L == 0 and RESOR < 1e-12)
	{
		cout<<EqnName<<": No need for inner iteration "<<endl;
		break;
	}
*/	
	else if(L == iterations - 1)
	{
		cout<<EqnName<<":SIP Solver Max iteration reached. Inner Iteration over "<<endl;

	}

	} // end iteration loop

	cout<<endl;
	forAll(phitemp)
	{
		phitemp[i][j][k].value = phi[i][j][k].value;
	}

	return phitemp;

} 

//------------------------------------ END OF SIP SOLVER FOR VELOCITY ----------------------------------------//




//---------------------------------------- SIP SOLVER FOR PRESSURE -------------------------------------------//

Fields::vec3dField Equation::sipSolverPressure(Fields::vec3dField& phi, Solution& sol, unsigned int& iterations,
						double& tolerance)
{
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), 
					Fields::vec1dField(phi[0][0].size())));

	// coefficients of upper and lower triangular matrix
	
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				LB[i][j][k].value = -AB[i][j][k].value
							/( 1.0 + sol.alfa*(UN[i][j][k-1].value 
								+ UE[i][j][k-1].value));
				LW[i][j][k].value = -AW[i][j][k].value
							/( 1.0 + sol.alfa*(UN[i-1][j][k].value 
								+ UT[i-1][j][k].value));
				LS[i][j][k].value = -AS[i][j][k].value
							/( 1.0 + sol.alfa*(UE[i][j-1][k].value 
								+ UT[i][j-1][k].value));

				double H1 = sol.alfa*(LB[i][j][k].value*UN[i][j][k-1].value
							+ LW[i][j][k].value*UN[i-1][j][k].value);
				double H2 = sol.alfa*(LB[i][j][k].value*UE[i][j][k-1].value
							+ LS[i][j][k].value*UE[i][j-1][k].value);
				double H3 = sol.alfa*(LW[i][j][k].value*UT[i-1][j][k].value
							+ LS[i][j][k].value*UT[i][j-1][k].value);

			
				LPR[i][j][k].value =1.0/( AP[i][j][k].value + H1 + H2 + H3 
				- LB[i][j][k].value*UT[i][j][k-1].value
				- LW[i][j][k].value*UE[i-1][j][k].value - LS[i][j][k].value*UN[i][j-1][k].value
				+ 1e-30);


				UN[i][j][k].value = (-AN[i][j][k].value - H1)*LPR[i][j][k].value;
				UE[i][j][k].value = (-AE[i][j][k].value - H2)*LPR[i][j][k].value;
				UT[i][j][k].value = (-AT[i][j][k].value - H3)*LPR[i][j][k].value;
			}

		}

	}
	

	//calculating and over-riding the residual
	for(unsigned int L = 0; L<iterations; L++)
	{
	Residual = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				RES[i][j][k].value = SP[i][j][k].value
				- (AP[i][j][k].value*phi[i][j][k].value
				- AE[i][j][k].value*phi[i+1][j][k].value 
				- AW[i][j][k].value*phi[i-1][j][k].value
				- AS[i][j][k].value*phi[i][j-1][k].value 
				- AN[i][j][k].value*phi[i][j+1][k].value
				- AT[i][j][k].value*phi[i][j][k+1].value 
				- AB[i][j][k].value*phi[i][j][k-1].value);

				Residual += std::pow(RES[i][j][k].value, 2.0);

				AUX[i][j][k].value = (RES[i][j][k].value 
						- LB[i][j][k].value*AUX[i][j][k-1].value
						- LW[i][j][k].value*AUX[i-1][j][k].value 
						- LS[i][j][k].value*AUX[i][j-1][k].value)*LPR[i][j][k].value; 
						// U*delta
				
				//LUdelta = RES
				//U*delta = Linv* RES
			}

		}

	}

	R2 = std::sqrt(Residual);
	double small = 1e-30;
	if(L==0)
	{
		RESOR = R2;
	}

	R2norm = R2/(RESOR + small);
	

	//Back substitution and correction
	for(unsigned int k = phi[0][0].size() - 2; k >= 1; --k)
	{
		for(unsigned int i = phi.size() - 2; i >= 1; --i)
		{
			for(unsigned int j = phi[0].size() - 2; j >= 1; --j)
			{
				DELTA[i][j][k].value = AUX[i][j][k].value 
						- UN[i][j][k].value*DELTA[i][j+1][k].value
						- UE[i][j][k].value*DELTA[i+1][j][k].value
						- UT[i][j][k].value*DELTA[i][j][k+1].value; 

				phi[i][j][k].value = phi[i][j][k].value + DELTA[i][j][k].value;
			}


		}

	}
	
	// check convergence
	cout<<EqnName<< " Inner iteration of SIP Solver "<< L+1 
				<< " and R2 --> " <<R2<< " R2norm "<< R2norm <<endl;
	if(R2norm <= SOR)
	{
		cout<<EqnName<<": Inner iteration of SIP Solver converged at iteration "<< L+1<<endl;
		break;
	}
	
/*
	if(R2norm <= tolerance)
	{
		cout<<EqnName<<": Inner iteration of SIP Solver converged at iteration "<< L+1<<endl;
		break;
	}
	
	if(L == 0 and RESOR < 1e-10)
	{
		cout<<EqnName<<": No need for inner iteration "<<endl;
		break;
	}

*/
	else if(L == iterations - 1)
	{
		cout<<EqnName<<":SIP Solver Max iteration reached. Inner iteration over "<<endl;

	}


	} // end iteration loop

	cout<<endl;
	forAll(phitemp)
	{
		phitemp[i][j][k].value = phi[i][j][k].value;
	}

	return phitemp;

} 

//------------------------------- END OF SIP SOLVER FOR PRESSURE ---------------------------------------------//




//-------------------------------- CGSTAB SOLVER FOR PRESSURE ------------------------------------------------//

Fields::vec3dField Equation::CGSTAB(Fields::vec3dField& phi, Solution& sol, 
				unsigned int& iterations, double& tolerance)
{
	
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), 
					Fields::vec1dField(phi[0][0].size())));

	// CALCULATION OF INITIAL RESIDUAL VECTOR

	double RES0 = 0.0;
	
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				RES[i][j][k].value = SP[i][j][k].value
				- (AP[i][j][k].value*phi[i][j][k].value
				- AE[i][j][k].value*phi[i+1][j][k].value 
				- AW[i][j][k].value*phi[i-1][j][k].value
				- AS[i][j][k].value*phi[i][j-1][k].value 
				- AN[i][j][k].value*phi[i][j+1][k].value
				- AT[i][j][k].value*phi[i][j][k+1].value 
				- AB[i][j][k].value*phi[i][j][k-1].value);

				RES0 += std::pow(RES[i][j][k].value, 2.0);
				RESO[i][j][k].value = RES[i][j][k].value;
			}

		}

	}
	RES0 = std::sqrt(RES0);
	RESOR = RES0;
	cout<<EqnName<< " Initial Residual of CGSTAB Solver " << " and R2 --> " <<RESOR <<endl;
	//CALCULATE ELEMENTS OF PRECONDITIONING MATRIX DIAGONAL

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				D[i][j][k].value = 1.0/(AP[i][j][k].value 
							- AW[i][j][k].value*D[i-1][j][k].value
							*AE[i-1][j][k].value
							- AS[i][j][k].value*D[i][j-1][k].value
							*AN[i][j-1][k].value
							- AB[i][j][k].value*D[i][j][k-1].value
							*AT[i][j][k-1].value);

			}

		}

	}

	// START OF ITERATIONS

	for(unsigned int L = 0; L<iterations; L++)
	{

	double BET = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				BET = BET + RES[i][j][k].value*RESO[i][j][k].value;
			}

		}

	}
	double OM = BET*GAM/(ALF*BETO + 1e-30);
	BETO=BET;

	// CALCULATE PK

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				PK[i][j][k].value = RES[i][j][k].value + OM*(PK[i][j][k].value 
									- ALF*UK[i][j][k].value);
			}

		}

	}

	//SOLVE (M ZK = PK) - FORWARD SUBSTITUTION

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = (PK[i][j][k].value 
							+ AW[i][j][k].value*ZK[i-1][j][k].value
							+ AS[i][j][k].value*ZK[i][j-1][k].value
							+ AB[i][j][k].value*ZK[i][j][k-1].value)
							*D[i][j][k].value;
			}

		}

	}

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = ZK[i][j][k].value/(D[i][j][k].value + 1e-30);
			}

		}

	}

	//BACKWARD SUBSTUITION

	for(unsigned int k = phi[0][0].size() - 2; k >= 1; --k)
	{
		for(unsigned int i = phi.size() - 2; i >= 1; --i)
		{
			for(unsigned int j = phi[0].size() - 2; j >= 1; --j)
			{
				ZK[i][j][k].value = (ZK[i][j][k].value 
							+ AE[i][j][k].value*ZK[i+1][j][k].value
							+ AN[i][j][k].value*ZK[i][j+1][k].value
							+ AT[i][j][k].value*ZK[i][j][k+1].value)
							*D[i][j][k].value;
			
			}


		}

	}
	
	//CALCULATE UK = A.PK

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				UK[i][j][k].value = AP[i][j][k].value*ZK[i][j][k].value
							- AW[i][j][k].value*ZK[i-1][j][k].value
							- AE[i][j][k].value*ZK[i+1][j][k].value
							- AN[i][j][k].value*ZK[i][j+1][k].value
							- AS[i][j][k].value*ZK[i][j-1][k].value
							- AT[i][j][k].value*ZK[i][j][k+1].value
							- AB[i][j][k].value*ZK[i][j][k-1].value;
			}

		}

	}

	//CALCULATE SCALAR PRODUCT UK.RESO AND GAMMA

	double UKRESO = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				UKRESO = UKRESO + UK[i][j][k].value*RESO[i][j][k].value;
			}

		}

	}
	GAM = BET/UKRESO;

	//UPDATE PHI AND CALCULATE W (OVERWRITE RES - IT IS RES-UPDATE)

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				phi[i][j][k].value = phi[i][j][k].value + GAM*ZK[i][j][k].value;
				RES[i][j][k].value = RES[i][j][k].value - GAM*UK[i][j][k].value;
			}

		}

	}

	//SOLVE (M Y = W) ; Y OVERWRITES ZK; FORWARD SUBSTITUTION

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = (RES[i][j][k].value 
							+ AW[i][j][k].value*ZK[i-1][j][k].value
							+ AS[i][j][k].value*ZK[i][j-1][k].value
							+ AB[i][j][k].value*ZK[i][j][k-1].value)
							*D[i][j][k].value;
			}

		}

	}

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = ZK[i][j][k].value/(D[i][j][k].value + 1e-30);
			}

		}

	}

	//BACKWARD SUBSTUITION

	for(unsigned int k = phi[0][0].size() - 2; k >= 1; --k)
	{
		for(unsigned int i = phi.size() - 2; i >= 1; --i)
		{
			for(unsigned int j = phi[0].size() - 2; j >= 1; --j)
			{
				ZK[i][j][k].value = (ZK[i][j][k].value 
							+ AE[i][j][k].value*ZK[i+1][j][k].value
							+ AN[i][j][k].value*ZK[i][j+1][k].value
							+ AT[i][j][k].value*ZK[i][j][k+1].value)
							*D[i][j][k].value;
			
			}


		}

	}

	//CALCULATE V = A.Y (VK = A.ZK)

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				VK[i][j][k].value = AP[i][j][k].value*ZK[i][j][k].value
							- AW[i][j][k].value*ZK[i-1][j][k].value
							- AE[i][j][k].value*ZK[i+1][j][k].value
							- AN[i][j][k].value*ZK[i][j+1][k].value
							- AS[i][j][k].value*ZK[i][j-1][k].value
							- AT[i][j][k].value*ZK[i][j][k+1].value
							- AB[i][j][k].value*ZK[i][j][k-1].value;
			}

		}

	}

	//CALCULATE ALPHA (ALF)

	double VRES = 0.0;
	double VV = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				VRES = VRES + VK[i][j][k].value*RES[i][j][k].value;
				VV = VV + VK[i][j][k].value*VK[i][j][k].value;
			}

		}

	}
	ALF = VRES/(VV + 1e-30);

	//UPDATE VARIABLE PHI AND RESIDUAL (RES) VECTORS

	double RESL = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				phi[i][j][k].value = phi[i][j][k].value + ALF*ZK[i][j][k].value;
				RES[i][j][k].value = RES[i][j][k].value - ALF*VK[i][j][k].value;
				RESL = RESL + std::pow(RES[i][j][k].value, 2.0);
			}

		}

	}
/*
	// TESTING ONLY

	cout<<" ITERATION NO: "<<L + 1<<endl;
	cout<<EqnName<<"  - RES " <<endl;
	testmat.print3dmat(RES);
	cout<<" END RES "<<endl;
	// TESTING END
*/
	R2 = std::sqrt(RESL);
	
	
	// CHECK CONVERGENCE

	R2norm = R2/(RESOR + 1e-30);
	
	cout<<EqnName<< " Inner iteration of CGSTAB Solver "<< L+1 
					<< " and R2 --> " <<R2<< " R2norm "<< R2norm <<endl;


	if(L == iterations - 1)
	{
		cout<<EqnName<<":CGSTAB Solver Max iteration reached. Inner iteration over "<<endl;
	}

	
	if(R2norm <= SOR)
	{
		cout<<EqnName<<": Inner iteration of CGSTAB Solver converged at iteration "<< L+1<<endl;
		break;
	}


	} // end iteration loop

	cout<<endl;
	forAll(phitemp)
	{
		phitemp[i][j][k].value = phi[i][j][k].value;
	}

	return phitemp;
}

//--------------------------------- END OF CGSTAB SOLVER FOR PRESSURE ----------------------------------------//




//------------------------------------CGSTAB SOLVER FOR VELOCITY ---------------------------------------------//

Fields::vec3dField Equation::CGSTABvel(Fields::vec3dField& phi, Fields::vec3dField& phiold, Solution& sol, 
						unsigned int& iterations, double& tolerance)
{
	
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), 
						Fields::vec1dField(phi[0][0].size())));

	// CALCULATION OF INITIAL RESIDUAL VECTOR

	double RES0 = 0.0;
	
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				RES[i][j][k].value = sourceRelaxed[i][j][k].value
				+ APNot[i][j][k].value*phiold[i][j][k].value
				- (AP[i][j][k].value*phi[i][j][k].value
				- AE[i][j][k].value*phi[i+1][j][k].value 
				- AW[i][j][k].value*phi[i-1][j][k].value
				- AS[i][j][k].value*phi[i][j-1][k].value 
				- AN[i][j][k].value*phi[i][j+1][k].value
				- AT[i][j][k].value*phi[i][j][k+1].value 
				- AB[i][j][k].value*phi[i][j][k-1].value);
				
				assert(isnan(RES[i][j][k].value) != 1);
				RES0 += std::pow(RES[i][j][k].value, 2.0);
				RESO[i][j][k].value = RES[i][j][k].value;
			}

		}

	}
	RES0 = std::sqrt(RES0);
	RESOR = RES0;
	cout<<EqnName<< " Initial Residual of CGSTAB Solver "<< " and R2 --> " <<RESOR <<endl;

	//CALCULATE ELEMENTS OF PRECONDITIONING MATRIX DIAGONAL

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				D[i][j][k].value = 1.0/(AP[i][j][k].value 
							- AW[i][j][k].value*D[i-1][j][k].value
							*AE[i-1][j][k].value
							- AS[i][j][k].value*D[i][j-1][k].value
							*AN[i][j-1][k].value
							- AB[i][j][k].value*D[i][j][k-1].value
							*AT[i][j][k-1].value);

			}

		}

	}


	// START OF ITERATIONS

	for(unsigned int L = 0; L<iterations; L++)
	{

	double BET = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				BET = BET + RES[i][j][k].value*RESO[i][j][k].value;
			}

		}

	}
	double OM = BET*GAM/(ALF*BETO + 1e-30);
	BETO=BET;

	// CALCULATE PK

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				PK[i][j][k].value = RES[i][j][k].value + OM*(PK[i][j][k].value 
									- ALF*UK[i][j][k].value);
			}

		}

	}

	//SOLVE (M ZK = PK) - FORWARD SUBSTITUTION

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = (PK[i][j][k].value 
							+ AW[i][j][k].value*ZK[i-1][j][k].value
							+ AS[i][j][k].value*ZK[i][j-1][k].value
							+ AB[i][j][k].value*ZK[i][j][k-1].value)
							*D[i][j][k].value;
			}

		}

	}

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = ZK[i][j][k].value/(D[i][j][k].value + 1e-30);
			}

		}

	}

	//BACKWARD SUBSTUITION

	for(unsigned int k = phi[0][0].size() - 2; k >= 1; --k)
	{
		for(unsigned int i = phi.size() - 2; i >= 1; --i)
		{
			for(unsigned int j = phi[0].size() - 2; j >= 1; --j)
			{
				ZK[i][j][k].value = (ZK[i][j][k].value 
							+ AE[i][j][k].value*ZK[i+1][j][k].value
							+ AN[i][j][k].value*ZK[i][j+1][k].value
							+ AT[i][j][k].value*ZK[i][j][k+1].value)
							*D[i][j][k].value;
			
			}


		}

	}
	
	//CALCULATE UK = A.PK

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				UK[i][j][k].value = AP[i][j][k].value*ZK[i][j][k].value
							- AW[i][j][k].value*ZK[i-1][j][k].value
							- AE[i][j][k].value*ZK[i+1][j][k].value
							- AN[i][j][k].value*ZK[i][j+1][k].value
							- AS[i][j][k].value*ZK[i][j-1][k].value
							- AT[i][j][k].value*ZK[i][j][k+1].value
							- AB[i][j][k].value*ZK[i][j][k-1].value;
			}

		}

	}

	//CALCULATE SCALAR PRODUCT UK.RESO AND GAMMA

	double UKRESO = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				UKRESO = UKRESO + UK[i][j][k].value*RESO[i][j][k].value;
			}

		}

	}
	GAM = BET/UKRESO;

	//UPDATE PHI AND CALCULATE W (OVERWRITE RES - IT IS RES-UPDATE)

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				phi[i][j][k].value = phi[i][j][k].value + GAM*ZK[i][j][k].value;
				RES[i][j][k].value = RES[i][j][k].value - GAM*UK[i][j][k].value;
			}

		}

	}

	//SOLVE (M Y = W) ; Y OVERWRITES ZK; FORWARD SUBSTITUTION

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = (RES[i][j][k].value 
							+ AW[i][j][k].value*ZK[i-1][j][k].value
							+ AS[i][j][k].value*ZK[i][j-1][k].value
							+ AB[i][j][k].value*ZK[i][j][k-1].value)
							*D[i][j][k].value;
			}

		}

	}

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = ZK[i][j][k].value/(D[i][j][k].value + 1e-30);
			}

		}

	}

	//BACKWARD SUBSTUITION

	for(unsigned int k = phi[0][0].size() - 2; k >= 1; --k)
	{
		for(unsigned int i = phi.size() - 2; i >= 1; --i)
		{
			for(unsigned int j = phi[0].size() - 2; j >= 1; --j)
			{
				ZK[i][j][k].value = (ZK[i][j][k].value + 
							AE[i][j][k].value*ZK[i+1][j][k].value
							+ AN[i][j][k].value*ZK[i][j+1][k].value
							+ AT[i][j][k].value*ZK[i][j][k+1].value)
							*D[i][j][k].value;
			
			}


		}

	}

	//CALCULATE V = A.Y (VK = A.ZK)

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				VK[i][j][k].value = AP[i][j][k].value*ZK[i][j][k].value
							- AW[i][j][k].value*ZK[i-1][j][k].value
							- AE[i][j][k].value*ZK[i+1][j][k].value
							- AN[i][j][k].value*ZK[i][j+1][k].value
							- AS[i][j][k].value*ZK[i][j-1][k].value
							- AT[i][j][k].value*ZK[i][j][k+1].value
							- AB[i][j][k].value*ZK[i][j][k-1].value;
			}

		}

	}

	//CALCULATE ALPHA (ALF)

	double VRES = 0.0;
	double VV = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				VRES = VRES + VK[i][j][k].value*RES[i][j][k].value;
				VV = VV + VK[i][j][k].value*VK[i][j][k].value;
			}

		}

	}
	ALF = VRES/(VV + 1e-30);

	//UPDATE VARIABLE PHI AND RESIDUAL (RES) VECTORS

	double RESL = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				phi[i][j][k].value = phi[i][j][k].value + ALF*ZK[i][j][k].value;
				RES[i][j][k].value = RES[i][j][k].value - ALF*VK[i][j][k].value;
				RESL = RESL + std::pow(RES[i][j][k].value, 2.0);
				assert(isnan(RES[i][j][k].value) != 1);
			}

		}

	}

/*	// TESTING ONLY

	cout<<" ITERATION NO: "<<iterations + 1<<endl;
	cout<<EqnName<<"  - RES " <<endl;
	testmat.print3dmat(RES);
	cout<<" END RES "<<endl;
	// TESTING END
*/

	R2 = std::sqrt(RESL);
	
	// CHECK CONVERGENCE

	R2norm = R2/(RESOR + 1e-30);
	
	cout<<EqnName<< " Inner iteration of CGSTAB Solver "<< L+1
				<< " and R2 --> " <<R2<< " R2norm "<< R2norm <<endl;


	if(L == iterations - 1)
	{
		cout<<EqnName<<":CGSTAB Solver Max iteration reached. Inner ieration over  "<<endl;

	}

	if(R2norm <= SOR)
	{
		cout<<EqnName<<": Inner iteration of CGSTAB Solver converged at iteration "<< L+1<<endl;
		break;
	}

	} // end iteration loop

	cout<<endl;
	forAll(phitemp)
	{
		phitemp[i][j][k].value = phi[i][j][k].value;
	}

	return phitemp;
}

//---------------------------------- END OF CGSTAB FOR VELOCITY ----------------------------------------------//





//-----------------------------------------ICCG FOR PRESSURE--------------------------------------------------//

Fields::vec3dField Equation::ICCG(Fields::vec3dField& phi, Solution& sol, 
unsigned int& iterations, double& tolerance)
{
	Fields::vec3dField phitemp (phi.size(), Fields::vec2dField(phi[0].size(), 
					Fields::vec1dField(phi[0][0].size())));
	
	double RES0 = 0.0;
	
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				RES[i][j][k].value = SP[i][j][k].value
				- (AP[i][j][k].value*phi[i][j][k].value
				- AE[i][j][k].value*phi[i+1][j][k].value 
				- AW[i][j][k].value*phi[i-1][j][k].value
				- AS[i][j][k].value*phi[i][j-1][k].value 
				- AN[i][j][k].value*phi[i][j+1][k].value
				- AT[i][j][k].value*phi[i][j][k+1].value 
				- AB[i][j][k].value*phi[i][j][k-1].value);

				RES0 += std::pow(RES[i][j][k].value, 2.0);
			}

		}

	}
	RES0 = std::sqrt(RES0);
	RESOR = RES0;

	//CALCULATE ELEMENTS OF PRECONDITIONING MATRIX DIAGONAL

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				D[i][j][k].value = 1.0/(AP[i][j][k].value 
							- std::pow(AW[i][j][k].value, 2.0)*D[i-1][j][k].value
							- std::pow(AS[i][j][k].value, 2.0)*D[i][j-1][k].value
							- std::pow(AB[i][j][k].value, 2.0)*D[i][j][k-1].value);

			}

		}

	}
	double S0 = 1e20;

	// START OF ITERATIONS

	for(unsigned int L = 0; L<iterations; L++)
	{

	//SOLVE FOR ZK - FORWARD SUBSTITUTION

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = (RES[i][j][k].value + AW[i][j][k].value*ZK[i-1][j][k].value
							+ AS[i][j][k].value*ZK[i][j-1][k].value
							+ AB[i][j][k].value*ZK[i][j][k-1].value)
							*D[i][j][k].value;
			}

		}

	}

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = ZK[i][j][k].value/(D[i][j][k].value + 1e-30);
			}

		}

	}

	//BACKWARD SUBSTUITION; CALCULATE SCALAR PRODUCT SK

	double SK = 0.0;
	for(unsigned int k = phi[0][0].size() - 2; k >= 1; --k)
	{
		for(unsigned int i = phi.size() - 2; i >= 1; --i)
		{
			for(unsigned int j = phi[0].size() - 2; j >= 1; --j)
			{
				ZK[i][j][k].value = (ZK[i][j][k].value
							+ AE[i][j][k].value*ZK[i+1][j][k].value
							+ AN[i][j][k].value*ZK[i][j+1][k].value
							+ AT[i][j][k].value*ZK[i][j][k+1].value)
							*D[i][j][k].value;
			
				SK = SK + RES[i][j][k].value*ZK[i][j][k].value;
			}


		}

	}


	//CALCULATE BETA

	double BET = SK/S0;


	//CALCULATE NEW SEARCH VECTOR PK

	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				PK[i][j][k].value = ZK[i][j][k].value + BET*PK[i][j][k].value;
			}

		}

	}

	//CALCULATE SCALAR PRODUCT (PK.A PK) AND ALPHA (OVERWRITE ZK)

	double PKAPK = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				ZK[i][j][k].value = AP[i][j][k].value*PK[i][j][k].value 
						- AE[i][j][k].value*PK[i+1][j][k].value
						- AW[i][j][k].value*PK[i-1][j][k].value
						- AN[i][j][k].value*PK[i][j+1][k].value
						- AS[i][j][k].value*PK[i][j-1][k].value
						- AT[i][j][k].value*PK[i][j][k+1].value
						- AB[i][j][k].value*PK[i][j][k-1].value;
				PKAPK = PKAPK + PK[i][j][k].value*ZK[i][j][k].value;
			}

		}

	}
	ALF = SK/PKAPK;

	//CALCULATE VARIABLE CORRECTION, NEW RESIDUAL VECTOR, AND NORM

	double RESL = 0.0;
	for(unsigned int k = 1; k< phi[0][0].size() -1; k++)
	{
		for(unsigned int i = 1; i<phi.size() -1; i++)
		{
			for(unsigned int j = 1; j<phi[0].size() - 1; j++)
			{
				phi[i][j][k].value = phi[i][j][k].value + ALF*PK[i][j][k].value;
				RES[i][j][k].value = RES[i][j][k].value - ALF*ZK[i][j][k].value;
				RESL = RESL + std::pow(RES[i][j][k].value, 2.0);
			}

		}

	}
	RESL = std::sqrt(RESL);
	R2 = RESL;
	S0 = SK;
	// CHECK CONVERGENCE

	R2norm = R2/(RESOR + 1e-30);
	
	cout<<EqnName<< " Inner iteration of ICCG Solver "<< L+1 << " and R2 --> " <<R2<< " R2norm "
						<< R2norm <<endl;
	if(L == 0 and RESOR < 1e-9)
	{
		cout<<EqnName<<": No need for inner iteration "<<endl;
		break;
	}

	if(L == iterations - 1)
	{
		cout<<EqnName<<":ICCG Solver Max iteration reached. NOT CONVERGED. EXITING PROGRAM "<<endl;
		exit(0);
	}

	if(R2 > 1e-9)
	{
		if(R2norm <= tolerance and R2 < 1e-6)
		{
			cout<<EqnName<<": Inner iteration of ICCG Solver converged at iteration "<< L+1<<endl;
			break;
		}
	}

	else
	{
		cout<<EqnName<<": Inner iteration of ICCG Solver converged at iteration "<< L+1<<endl;
		break;

	}


	} // end iteration loop

	cout<<endl;
	forAll(phitemp)
	{
		if(std::abs(phi[i][j][k].value)<1e-8)
		{
			phi[i][j][k].value = 0.0;
		}

		phitemp[i][j][k].value = phi[i][j][k].value;
	}

	return phitemp;

}

//--------------------------------END OF ICCG FOR PRESSURE----------------------------------------------------//






//------------------------- RHIE - CHOW MOMENTUM INTERPOLATION -----------------------------------------------//

Fields::vec3dField Equation::momentumInterpolation(Fields::vec3dField& vel, Fields::vec3dField& pr, 
		Fields::vec3dField& velprev, Fields::vec3dField& velWallold, Fields::vec3dField& velold,
		Fields::vec3dField& UvelWallprev, Fields::vec3dField& VvelWallprev,
		Fields::vec3dField& WvelWallprev, Solution& sol_, int& direction)
{
		
	if(direction == 1) // U vel interpolation
	{
		Fields::vec3dField wallvel (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
		
		forAllInternal(wallvel)
		{
			double f = vel[i][j][k].FXE;

			double HP = (AE[i][j][k].value*vel[i+1][j][k].value 
				+ AW[i][j][k].value*vel[i-1][j][k].value
				+ AS[i][j][k].value*vel[i][j-1][k].value 
				+ AN[i][j][k].value*vel[i][j+1][k].value
				+ AT[i][j][k].value*vel[i][j][k+1].value 
				+ AB[i][j][k].value*vel[i][j][k-1].value
				+ (1 - URF)*velprev[i][j][k].value*AP[i][j][k].value)*rAP[i][j][k].value;

			double HE = (AE[i+1][j][k].value*vel[i+2][j][k].value 
				+ AW[i+1][j][k].value*vel[i][j][k].value
				+ AS[i+1][j][k].value*vel[i+1][j-1][k].value 
				+ AN[i+1][j][k].value*vel[i+1][j+1][k].value
				+ AT[i+1][j][k].value*vel[i+1][j][k+1].value 
				+ AB[i+1][j][k].value*vel[i+1][j][k-1].value
				+ (1 - URF)*velprev[i+1][j][k].value*AP[i+1][j][k].value)*rAP[i+1][j][k].value;
			

			double H_ = f*HE + (1-f)*HP;

			double AP_ = AP[i+1][j][k].value*AP[i][j][k].value
							/(f*AP[i][j][k].value + (1-f)*AP[i+1][j][k].value);
				
			double Press_BodyF_Time = 
				vel[i][j][k].volume*(pr[i][j][k].value - pr[i+1][j][k].value)
							/vel[i][j][k].DXPtoE
				+ sol_.density*velWallold[i][j][k].value*vel[i][j][k].volume/sol_.dt
				- 2*sol_.density*WvelWallprev[i][j][k].value*sol_.omega*vel[i][j][k].volume
				+ sol_.density*sol_.omega*sol_.omega*vel[i][j][k].Se
				*(vel[i+1][j][k].XC*vel[i+1][j][k].XC - vel[i][j][k].XC*vel[i][j][k].XC)/2;


			double wallvelrelaxed = (1-URF)*(UvelWallprev[i][j][k].value - 
							(f*velprev[i+1][j][k].value + (1-f)*velprev[i][j][k].value));

			double wallveloldrelaxed = sol_.density*vel[i][j][k].volume/(AP_*sol_.dt)
								*(velWallold[i][j][k].value - 
							(f*velold[i+1][j][k].value + (1-f)*velold[i][j][k].value));


			wallvel[i][j][k].value = H_ + Press_BodyF_Time/AP_ + wallvelrelaxed 
						+ wallveloldrelaxed;


		}

		forWestBoundary(wallvel)
		{
			wallvel[i][j][k].value = 0.0;

		}	

		forEastBoundary(wallvel)
		{
			wallvel[i][j][k].value = 0.0;

		}
		return wallvel;
	}

	else if (direction == 2) // V vel interpolation
	{
		Fields::vec3dField wallvel (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));

		forAllInternal(wallvel)
		{
			double f = vel[i][j][k].FYN;

			double HP = (AE[i][j][k].value*vel[i+1][j][k].value 
				+ AW[i][j][k].value*vel[i-1][j][k].value
				+ AS[i][j][k].value*vel[i][j-1][k].value 
				+ AN[i][j][k].value*vel[i][j+1][k].value
				+ AT[i][j][k].value*vel[i][j][k+1].value 
				+ AB[i][j][k].value*vel[i][j][k-1].value
				+ (1 - URF)*velprev[i][j][k].value*AP[i][j][k].value)*rAP[i][j][k].value;

			double HN = (AE[i][j+1][k].value*vel[i+1][j+1][k].value 
				+ AW[i][j+1][k].value*vel[i-1][j+1][k].value
				+ AS[i][j+1][k].value*vel[i][j][k].value 
				+ AN[i][j+1][k].value*vel[i][j+2][k].value
				+ AT[i][j+1][k].value*vel[i][j+1][k+1].value 
				+ AB[i][j+1][k].value*vel[i][j+1][k-1].value
				+ (1 - URF)*velprev[i][j+1][k].value*AP[i][j+1][k].value)*rAP[i][j+1][k].value;

			double H_ = f*HN + (1-f)*HP;

			double AP_ = AP[i][j+1][k].value*AP[i][j][k].value
					/(f*AP[i][j][k].value + (1-f)*AP[i][j+1][k].value);
				
			double Press_BodyF_Time = 
				vel[i][j][k].volume*(pr[i][j][k].value - pr[i][j+1][k].value)
							/vel[i][j][k].DYPtoN
				+ sol_.density*velWallold[i][j][k].value*vel[i][j][k].volume/sol_.dt;

			double wallvelrelaxed = (1-URF)*(VvelWallprev[i][j][k].value - 
							(f*velprev[i][j+1][k].value + (1-f)*velprev[i][j][k].value));

			double wallveloldrelaxed = sol_.density*vel[i][j][k].volume/(AP_*sol_.dt)
								*(velWallold[i][j][k].value - 
							(f*velold[i][j+1][k].value + (1-f)*velold[i][j][k].value));


			wallvel[i][j][k].value = H_ + Press_BodyF_Time/AP_ + wallvelrelaxed 
						+ wallveloldrelaxed;

		}

		
		forSouthBoundary(wallvel)
		{
			wallvel[i][j][k].value = 0.0;
		}
		
		forNorthBoundary(wallvel)
		{
		
			wallvel[i][j][k].value = 0.0;
		}
		
		return wallvel;

	}

	else if (direction == 3) // W vel interpolation
	{
		Fields::vec3dField wallvel (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));

	
		forAllInternal(wallvel)
		{
			double f = vel[i][j][k].FZT;

			double HP = (AE[i][j][k].value*vel[i+1][j][k].value 
				+ AW[i][j][k].value*vel[i-1][j][k].value
				+ AS[i][j][k].value*vel[i][j-1][k].value 
				+ AN[i][j][k].value*vel[i][j+1][k].value
				+ AT[i][j][k].value*vel[i][j][k+1].value 
				+ AB[i][j][k].value*vel[i][j][k-1].value
				+ (1 - URF)*velprev[i][j][k].value*AP[i][j][k].value)*rAP[i][j][k].value;


			double HT = (AE[i][j][k+1].value*vel[i+1][j][k+1].value 
				+ AW[i][j][k+1].value*vel[i-1][j][k+1].value
				+ AS[i][j][k+1].value*vel[i][j-1][k+1].value 
				+ AN[i][j][k+1].value*vel[i][j+1][k+1].value
				+ AT[i][j][k+1].value*vel[i][j][k+2].value 
				+ AB[i][j][k+1].value*vel[i][j][k].value
				+ (1 - URF)*velprev[i][j][k+1].value*AP[i][j][k+1].value)*rAP[i][j][k+1].value;

			double H_ = f*HT + (1-f)*HP;

			double AP_ = AP[i][j][k+1].value*AP[i][j][k].value
						/(f*AP[i][j][k].value + (1-f)*AP[i][j][k+1].value);
				
			double Press_BodyF_Time =
				vel[i][j][k].volume*(pr[i][j][k].value - pr[i][j][k+1].value)
						/vel[i][j][k].DZPtoT
				+ sol_.density*velWallold[i][j][k].value*vel[i][j][k].volume/sol_.dt
				+ 2*sol_.density*UvelWallprev[i][j][k].value*sol_.omega*vel[i][j][k].volume
				+ sol_.density*sol_.omega*sol_.omega*vel[i][j][k].Se
				*(vel[i][j][k+1].ZC*vel[i][j][k+1].ZC - vel[i][j][k].ZC*vel[i][j][k].ZC)/2;
			
			double wallvelrelaxed = (1-URF)*(WvelWallprev[i][j][k].value - 
							(f*velprev[i][j][k+1].value + (1-f)*velprev[i][j][k].value));

			double wallveloldrelaxed = sol_.density*vel[i][j][k].volume/(AP_*sol_.dt)
								*(velWallold[i][j][k].value - 
							(f*velold[i][j][k+1].value + (1-f)*velold[i][j][k].value));


			wallvel[i][j][k].value = H_ + Press_BodyF_Time/AP_ + wallvelrelaxed 
						+ wallveloldrelaxed;

		}

		forBottomBoundary(wallvel)
		{
			wallvel[i][j][k].value = vel[i][j][k].value;
		}
			


		forTopBoundary(wallvel)
		{
			wallvel[i][j][k].value = vel[i][j][k+1].value;
//			wallvel[i][j][k].dt = vel[i][j][k].St*rAP[i][j][k].value/2.0;
		}


		return wallvel;
	}
	
}

//---------------------------- END OF RHIE-CHOW MOMENTUM INTERPOLATION --------------------------------------//




//------------------------------------MASS FLUX CALCULATION -------------------------------------------------//


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

//------------------------------- END OF MASS FLUX CALCULATION -----------------------------------------------//





//---------------------------------- FACE VELOCITY CORRECTION ------------------------------------------------//

Fields::vec3dField Equation::faceVelocityCorrection(Fields::vec3dField& PCorr, Fields::vec3dField& UWall, 
						Fields::vec3dField& VWall, Fields::vec3dField& WWall, 
						int& direction)
{
	if(direction == 1) // Uwall vel correction
	{
		Fields::vec3dField wallvel (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

		forAllInternal(wallvel)
		{
			double f = PCorr[i][j][k].FXE;
			double AP_ = AP[i+1][j][k].value*AP[i][j][k].value
				/(f*AP[i][j][k].value + (1-f)*AP[i+1][j][k].value);

			wallvel[i][j][k].value = PCorr[i][j][k].Se
				*(PCorr[i][j][k].value - PCorr[i+1][j][k].value)/AP_;
		}

		return wallvel;
	}



	else if (direction == 2) // Vwall vel correction
	{
		Fields::vec3dField wallvel (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));


		forAllInternal(wallvel)
		{
			double f = PCorr[i][j][k].FYN;
			double AP_ = AP[i][j+1][k].value*AP[i][j][k].value
				/(f*AP[i][j][k].value + (1-f)*AP[i][j+1][k].value);

			wallvel[i][j][k].value = PCorr[i][j][k].Sn
				*(PCorr[i][j][k].value - PCorr[i][j+1][k].value)/AP_;

		}

		return wallvel;


	}



	else if (direction == 3) // Wwall vel correction
	{
		Fields::vec3dField wallvel (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));


		forAllInternal(wallvel)
		{
			double f = PCorr[i][j][k].FZT;
			double AP_ = AP[i][j][k+1].value*AP[i][j][k].value
				/(f*AP[i][j][k].value + (1-f)*AP[i][j][k+1].value);

			wallvel[i][j][k].value = PCorr[i][j][k].St
				*(PCorr[i][j][k].value - PCorr[i][j][k+1].value)/AP_;
	
		}
	/*	
		forTopBoundary(wallvel)
		{
			wallvel[i][j][k].value = WWall[i][j][k].dt
			*(PCorr[i][j][k].value - PCorr[i][j][k+1].value);
		}
	*/

		return wallvel;
	}

}

//--------------------------------- END OF FACE VELOCITY CORRECTION -----------------------------------------//




//-------------------------------------- MASSFLUX CORRECTION ------------------------------------------------//

Fields::vec3dField Equation::massFluxCorrection(Fields::vec3dField& PCorr, Fields::vec3dField& UWCorr, 
						Fields::vec3dField& VWCorr, Fields::vec3dField& WWCorr, 
						int& direction)
{
	if(direction == 1) // east direction mass flux correction
	{
		Fields::vec3dField massFlux (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
	
		forAllInternal(massFlux)
		{
			massFlux[i][j][k].value =PCorr[i][j][k].Se*PCorr[i][j][k].density*UWCorr[i][j][k].value;
		}
		
		return massFlux;
	}

	else if (direction == 2) // north direction mass flux
	{
		Fields::vec3dField massFlux (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
	
		forAllInternal(massFlux)
		{
			massFlux[i][j][k].value =PCorr[i][j][k].Sn*PCorr[i][j][k].density*VWCorr[i][j][k].value;
		}
		
		return massFlux;

	}

	else if (direction == 3) // top direction mass flux
	{
		Fields::vec3dField massFlux (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));
		
		forAllInternal(massFlux)
		{
			massFlux[i][j][k].value =PCorr[i][j][k].St*PCorr[i][j][k].density*WWCorr[i][j][k].value;
		}

	/*	
		forTopBoundary(massFlux)
		{
			massFlux[i][j][k].value =PCorr[i][j][k].St*PCorr[i][j][k].density*WWCorr[i][j][k].value;
		}
	*/

		return massFlux;

	}

}
//-------------------------------END OF MASS FLUX CORRECTION--------------------------------------------------//




//--------------------------------CELL VELOCITY CORRECTION----------------------------------------------------//

Fields::vec3dField Equation::cellVelocityCorrection(Fields::vec3dField& PCorr, int& direction)
{

	
	if(direction == 1) // U vel correction
	{
		Fields::vec3dField Ucell (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

		forAllInternal(Ucell)
		{
			double DX =  PCorr[i][j][k].X - PCorr[i-1][j][k].X;

			double PCE = PCorr[i+1][j][k].value*PCorr[i][j][k].FXE
					+ PCorr[i][j][k].value*PCorr[i][j][k].FXP;

			double PCW = PCorr[i][j][k].value*PCorr[i-1][j][k].FXE
					+ PCorr[i-1][j][k].value*PCorr[i-1][j][k].FXP;

			Ucell[i][j][k].value = rAP[i][j][k].value*PCorr[i][j][k].volume*(PCW - PCE)/DX;
		}

		return Ucell;
	}

	else if (direction == 2) // V vel correction
	{
		Fields::vec3dField Vcell (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));


		forAllInternal(Vcell)
		{
			double DY = PCorr[i][j][k].Y - PCorr[i][j-1][k].Y;

			double PCN = PCorr[i][j+1][k].value*PCorr[i][j][k].FYN
					+ PCorr[i][j][k].value*PCorr[i][j][k].FYP;
					
			double PCS = PCorr[i][j][k].value*PCorr[i][j-1][k].FYN
					+ PCorr[i][j-1][k].value*PCorr[i][j-1][k].FYP;

			Vcell[i][j][k].value = rAP[i][j][k].value*PCorr[i][j][k].volume*(PCS - PCN)/DY;

		}

		return Vcell;


	}

	else if (direction == 3) // W vel correction
	{
		Fields::vec3dField Wcell (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));


		forAllInternal(Wcell)
		{
			double DZ = PCorr[i][j][k].Z - PCorr[i][j][k-1].Z;

			double PCT = PCorr[i][j][k+1].value*PCorr[i][j][k].FZT
					+ PCorr[i][j][k].value*PCorr[i][j][k].FZP;
					
			double PCB = PCorr[i][j][k].value*PCorr[i][j][k-1].FZT
					+ PCorr[i][j][k-1].value*PCorr[i][j][k-1].FZP;

			Wcell[i][j][k].value = rAP[i][j][k].value*PCorr[i][j][k].volume*(PCB - PCT)/DZ;
	
		}

		return Wcell;
	}

}

//---------------------------- END OF CELL VELOCITY CORRECTION -----------------------------------------------//




Equation::~Equation()
{


}
