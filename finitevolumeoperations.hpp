#include "mesh.hpp"
#include "fields.hpp"
#include "solution.hpp"
#include "foralloperations.hpp"


namespace fvm
{
//-----------------------------------------MOMENTUM EQUATION------------------------------------------------//

//-------------------------------------CONVECTION-DIFFUSION TERM--------------------------------------------//

FiniteMatrix::finiteMat convDiffusiveTerm(Fields::vec3dField& vec, 
Fields::vec3dField& massFluxEast, Fields::vec3dField& massFluxNorth, Fields::vec3dField& massFluxTop)
{
	// recieves the massfluxes and velocity field ( only for distance, viscosity and surface area)
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), 
							vector<FiniteMatrix>(vec[0][0].size())));
	
	//towards east side
	forAllInternalUCVs(vec)
	{
		double Fe = massFluxEast[i][j][k].value;
		double De = (vec[i][j][k].visc*vec[i][j][k].Se)/vec[i][j][k].DXPtoE;

		APTemp[i][j][k].ae     = max(max(-Fe , (De - Fe/2.0)), 0.0);
		APTemp[i+1][j][k].aw   = max(max( Fe , (De + Fe/2.0)), 0.0);
	}


	//towards north side
	forAllInternalVCVs(vec)
	{
		double Fn = massFluxNorth[i][j][k].value;
		double Dn = (vec[i][j][k].visc*vec[i][j][k].Sn)/vec[i][j][k].DYPtoN;
	
		APTemp[i][j][k].an     =  max(max(-Fn , (Dn - Fn/2.0)), 0.0);
		APTemp[i][j+1][k].as   =  max(max( Fn , (Dn + Fn/2.0)), 0.0);

	
	}

	//towards top side
	forAllInternalWCVs(vec)
	{
		double Ft = massFluxTop[i][j][k].value;
		double Dt = (vec[i][j][k].visc*vec[i][j][k].St)/vec[i][j][k].DZPtoT;

		APTemp[i][j][k].at = max(max(-Ft , (Dt - Ft/2.0)), 0.0);
		APTemp[i][j][k+1].ab = max(max( Ft , (Dt + Ft/2.0)), 0.0);
	
	}

	//for cell near outlet boundary (top)

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = vec[i][j].size() - 2; k < vec[i][j].size() - 1; k++)
			{

				APTemp[i][j][k].at = 0.0;
				
			}
		}
	
	}

	
	//for cell near inlet boundary (bottom)
	
	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < 2; k++)
			{
				double Ft =  massFluxTop[i][j][k].value;
				double Dt =  (vec[i][j][k].visc*vec[i][j][k].St)/vec[i][j][k].DZPtoT;
				
				APTemp[i][j][k].at = max(max(-Ft, (-Ft/2.0 + 4.0*Dt/3.0)), 0.0);
				APTemp[i][j][k].ab = 0.0;

			}
		}
	
	}


	//cell near wall north

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = vec[i].size() - 2; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				double Fs = massFluxNorth[i][j-1][k].value;
				double Ds = (vec[i][j][k].visc*vec[i][j][k].Sn)/vec[i][j-1][k].DYPtoN;
				
				APTemp[i][j][k].as = max(max(Fs, Fs/2.0 + 4.0*Ds/3.0), 0.0);
				APTemp[i][j][k].an =  0.0;

			}
		}

	}

	// cell near Wall south

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<2; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				double Fn = massFluxNorth[i][j][k].value;
				double Dn = (vec[i][j][k].visc*vec[i][j][k].Sn)/vec[i][j][k].DYPtoN;
			
				APTemp[i][j][k].an = max(max(-Fn, -Fn/2.0 + 4.0*Dn/3.0), 0.0);
				APTemp[i][j][k].as = 0.0;

			}
		}
	}
	//cell near Wall east

	for(unsigned int i = vec.size() - 2; i<vec.size() - 1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				double Fw = massFluxEast[i-1][j][k].value;
				double Dw = (vec[i][j][k].visc*vec[i][j][k].Se)/vec[i-1][j][k].DXPtoE;

				APTemp[i][j][k].aw = max(max(Fw, (Fw/2.0 + 4.0*Dw/3.0)), 0.0);
				APTemp[i][j][k].ae = 0.0;

			}
		}
	}

	//cell near Wall west

	for(unsigned int i = 1; i < 2; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				double Fe = massFluxEast[i][j][k].value;
				double De = (vec[i][j][k].visc*vec[i][j][k].Se)/vec[i][j][k].DXPtoE;

				APTemp[i][j][k].aw = 0.0; 
				APTemp[i][j][k].ae = max(max(-Fe, (-Fe/2.0 + 4.0*De/3.0)), 0.0);

			}
		}

	}
	return APTemp;

}
//----------------------------------END CONVECTION-DIFFUSION TERM---------------------------------------------//





//--------------------------------PRESSURE SOURCE AND BOUNDARY SOURCES----------------------------------------//

FiniteMatrix::finiteMat pressureGrad(Fields::vec3dField& vec, Fields::vec3dField& Ucell, 
Fields::vec3dField& Vcell, Fields::vec3dField& Wcell, Fields::vec3dField& massFluxTop, int& direction_) 
// receives the pressure field
{
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), 
						vector<FiniteMatrix>(vec[0][0].size())));
	
	forAllInternal(vec)
	{
				
        	double pressureEastFace = (vec[i+1][j][k].value*vec[i][j][k].FXE) 
					+ (vec[i][j][k].value* (1.0-vec[i][j][k].FXE));
        	double pressureWestFace = (vec[i][j][k].value*vec[i-1][j][k].FXE) 
					+ (vec[i-1][j][k].value*(1.0-vec[i-1][j][k].FXE));


        	double pressureNorthFace = (vec[i][j+1][k].value*vec[i][j][k].FYN) 
					+ (vec[i][j][k].value*(1.0-vec[i][j][k].FYN));
        	double pressureSouthFace = (vec[i][j][k].value*vec[i][j-1][k].FYN) 
					+ (vec[i][j-1][k].value*(1.0-vec[i][j-1][k].FYN));

		double pressureTopFace    = (vec[i][j][k+1].value*vec[i][j][k].FZT)
					+ (vec[i][j][k].value*(1.0-vec[i][j][k].FZT));
        	double pressureBottomFace = (vec[i][j][k].value*vec[i][j][k-1].FZT) 
					+ (vec[i][j][k-1].value*(1.0-vec[i][j][k-1].FZT));

        	double pressureEastGrad  = (pressureEastFace - pressureWestFace);
       		double pressureNorthGrad = (pressureNorthFace- pressureSouthFace);
		double pressureTopGrad   = (pressureTopFace  - pressureBottomFace);

        	if(direction_ == 1)  //for east direction
        	{
            		APTemp[i][j][k].sp = -pressureEastGrad*vec[i][j][k].Se;

	       	}
		
        	else if(direction_ == 2) //for north direction
        	{
            		APTemp[i][j][k].sp = -pressureNorthGrad*vec[i][j][k].Sn;
		
        	}

		else if(direction_ == 3) // for top direction
		{
			APTemp[i][j][k].sp = -pressureTopGrad*vec[i][j][k].St;
			
		}

	}

        if(direction_ == 1)  //for east direction
        {     
		forBottomBoundary(APTemp)
		{
			APTemp[i][j][k+1].sp += (massFluxTop[i][j][k].value + 
						8.0*Ucell[i][j][k+1].visc*Ucell[i][j][k+1].St
						/(3.0*Ucell[i][j][k+1].DZPtoT))*Ucell[i][j][k].value;
				// source terms due to known inlet velocity
		}

       	}
		
        else if(direction_ == 2) //for north direction
        {
		forBottomBoundary(APTemp)
		{
			APTemp[i][j][k+1].sp += (massFluxTop[i][j][k].value + 
						8.0*Vcell[i][j][k+1].visc*Vcell[i][j][k+1].St
						/(3.0*Vcell[i][j][k+1].DZPtoT))*Vcell[i][j][k].value;
		}

        }

	else if(direction_ == 3) // for top direction
	{
	
		forBottomBoundary(APTemp)
		{
			APTemp[i][j][k+1].sp += (massFluxTop[i][j][k].value + 
						8.0*Wcell[i][j][k+1].visc*Wcell[i][j][k+1].St
						/(3.0*Wcell[i][j][k+1].DZPtoT))*Wcell[i][j][k].value;
		}
	}

	return APTemp;
}

//--------------------------------END PRESSURE SOURCE AND BOUNDARY SOURCES-------------------------------//




//--------------------------------------------FORCE SOURCE-----------------------------------------------//

// receives the velocity fields to calculate the coriolis force and centrifugal force

FiniteMatrix::finiteMat forceSource(Fields::vec3dField& vec, Solution& sol_, int& direction_) 
{
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), 
							vector<FiniteMatrix>(vec[0][0].size())));
	
	
	forAllInternal(vec)
	{
		if(direction_ == 1)  //for east direction (pass W velocity field )
        	{
            		APTemp[i][j][k].sf = -2*vec[i][j][k].density*vec[i][j][k].value*sol_.omega*vec[i][j][k].volume
					        + vec[i][j][k].density*sol_.omega*sol_.omega*vec[i][j][k].Se
					        *(vec[i][j][k].X*vec[i][j][k].X - vec[i-1][j][k].X*vec[i-1][j][k].X)/2;

       		}
		
        	else if(direction_ == 2) //for north direction
        	{
            		APTemp[i][j][k].sf = 0.0;
        	}

		else if(direction_ == 3) // for top direction (pass U velocity field)
		{
			APTemp[i][j][k].sf = 2*vec[i][j][k].density*vec[i][j][k].value*sol_.omega*vec[i][j][k].volume
					     + vec[i][j][k].density*sol_.omega*sol_.omega*vec[i][j][k].St
					     *(vec[i][j][k].Z*vec[i][j][k].Z - vec[i][j][k-1].Z*vec[i][j][k-1].Z)/2;
		}

		
        }
		
	return APTemp;
}

//--------------------------------------------END FORCE SOURCE----------------------------------------------//




//---------------------------------------------TIME COEFFICIENT---------------------------------------------//

FiniteMatrix::finiteMat timeCoeff(Fields::vec3dField& vec, Solution& sol_) 
{
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), 
							vector<FiniteMatrix>(vec[0][0].size())));
	
	
	forAllInternal(vec)
	{
		APTemp[i][j][k].apnot = sol_.density*vec[i][j][k].volume/sol_.dt;
		
        }
		
	return APTemp;
}

//----------------------------------------END TIME COEFFICIENT----------------------------------------------//





//----------------------------------PRESSURE CORRECTION EQN COEFFICIENT-------------------------------------//
// vec is pressure

FiniteMatrix::finiteMat pressureCorrectionCoeff(Fields::vec3dField& vec, Fields::vec3dField& massFluxEast, 
Fields::vec3dField& massFluxNorth, Fields::vec3dField& massFluxTop, FiniteMatrix::finiteMat& APU,
FiniteMatrix::finiteMat& APV, FiniteMatrix::finiteMat& APW)
{
	// recieves the massfluxes and velocity field ( only for distance, viscosity and surface area)
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), 
						vector<FiniteMatrix>(vec[0][0].size())));

	forAllInternalUCVs(vec)
	{
		double f = vec[i][j][k].FXE;
		double AP_ = APU[i+1][j][k].value*APU[i][j][k].value
				/(f*APU[i][j][k].value + (1-f)*APU[i+1][j][k].value);
	
		APTemp[i][j][k].ae = vec[i][j][k].density*vec[i][j][k].Se*vec[i][j][k].Se/AP_;
		APTemp[i+1][j][k].aw = vec[i+1][j][k].density*vec[i+1][j][k].Se*vec[i+1][j][k].Se/AP_;
	}

	forAllInternalVCVs(vec)
	{
		double f = vec[i][j][k].FYN;
		double AP_ = APV[i][j+1][k].value*APV[i][j][k].value
				/(f*APV[i][j][k].value + (1-f)*APV[i][j+1][k].value);

		APTemp[i][j][k].an = vec[i][j][k].density*vec[i][j][k].Sn*vec[i][j][k].Sn/AP_;
		APTemp[i][j+1][k].as = vec[i][j+1][k].density*vec[i][j+1][k].Sn*vec[i][j+1][k].Sn/AP_;
	}

	forAllInternalWCVs(vec)
	{

		double f = vec[i][j][k].FZT;
		double AP_ = APW[i][j][k+1].value*APW[i][j][k].value
				/(f*APW[i][j][k].value + (1-f)*APW[i][j][k+1].value);

		APTemp[i][j][k].at = vec[i][j][k].density*vec[i][j][k].St*vec[i][j][k].St/AP_;
		APTemp[i][j][k+1].ab = vec[i][j][k+1].density*vec[i][j][k+1].St*vec[i][j][k+1].St/AP_;
	}

	//for cell near outlet boundary (top)

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = vec[i][j].size() - 2; k < vec[i][j].size() - 1; k++)
			{
				APTemp[i][j][k].at = 0.0;
				APTemp[i][j][k].ab = APTemp[i][j][k].ab
			 - 0.5*vec[i][j][k].density*vec[i][j][k].St*vec[i][j][k].St/APW[i][j][k].value;
			}
		}
	
	}

	
	//for cell near inlet boundary (bottom)
	
	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < 2; k++)
			{

				APTemp[i][j][k].ab 	= 0.0;

			}
		}
	
	}


	//for cell near wall north

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = vec[i].size() - 2; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
	
				APTemp[i][j][k].an      =  0.0;

			}
		}

	}

	//for cell near Wall south

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<2; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
			
				APTemp[i][j][k].as   	=  0.0;

			}
		}
	}
	//for cell near Wall east

	for(unsigned int i = vec.size() - 2; i<vec.size() - 1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				
				APTemp[i][j][k].ae	= 0.0;
			}
		}
	}

	//for cell near Wall west

	for(unsigned int i = 1; i < 2; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				APTemp[i][j][k].aw	= 0.0;

			}
		}

	}

	
	return APTemp;

} 
//----------------------------------END OF PRESSURE CORRECTION COEFFICIENT------------------------------------//




//---------------------------------------------MASS IMBALANCE-------------------------------------------------//

FiniteMatrix::finiteMat pressureCorrectionSource(Fields::vec3dField& vec, Fields::vec3dField& massFluxEast, 
			Fields::vec3dField& massFluxNorth, Fields::vec3dField& massFluxTop)
{
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), 
							vector<FiniteMatrix>(vec[0][0].size())));
	//vec - pressure
	
	forAllInternal(vec)
	{
	
            	APTemp[i][j][k].sp = -massFluxEast[i][j][k].value 
				+ massFluxEast[i-1][j][k].value 
				- massFluxNorth[i][j][k].value 
				+ massFluxNorth[i][j-1][k].value 
				- massFluxTop[i][j][k].value 
				+ massFluxTop[i][j][k-1].value;
        }
		
	return APTemp;


}
//--------------------------------------END OF MASS IMBALANCE------------------------------------------------//

} 

//----------------------------------------END NAMESPACE FVM--------------------------------------------------//
