#include "mesh.hpp"
#include "fields.hpp"
#include "solution.hpp"
#include "foralloperations.hpp"


namespace fvm
{
//momentum equation 
//diffusion, convection, pressure source, time derivative, other sources

FiniteMatrix::finiteMat convDiffusiveTerm(Fields::vec3dField& vec, Fields::vec3dField& massFluxEast, Fields::vec3dField& massFluxNorth, Fields::vec3dField& massFluxTop)
{
	// recieves the massfluxes and velocity field ( only for distance, viscosity and surface area)
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), vector<FiniteMatrix>(vec[0][0].size())));
	
	//towards east side
	forAllInternalUCVs(vec)
	{
		double Fe = massFluxEast[i][j][k].value;
		double De = (vec[i][j][k].visc*vec[i][j][k].Se)/vec[i][j][k].DXPtoE;

		APTemp[i][j][k].ae     =  (De - Fe/2.0);//max(max(-Fe , (De - Fe/2.0)), 0.0);
		APTemp[i+1][j][k].aw   =  (De + Fe/2.0);//max(max( Fe , (De + Fe/2.0)), 0.0);
	}


	//towards north side
	forAllInternalVCVs(vec)
	{
		double Fn = massFluxNorth[i][j][k].value;
		double Dn = (vec[i][j][k].visc*vec[i][j][k].Sn)/vec[i][j][k].DYPtoN;
	
		APTemp[i][j][k].an     =  (Dn - Fn/2.0);//max(max(-Fn , (Dn - Fn/2.0)), 0.0);
		APTemp[i][j+1][k].as   =  (Dn + Fn/2.0);//max(max( Fn , (Dn + Fn/2.0)), 0.0);

	
	}

	//towards top side
	forAllInternalWCVs(vec)
	{
		double Ft = massFluxTop[i][j][k].value;
		double Dt = (vec[i][j][k].visc*vec[i][j][k].St)/vec[i][j][k].DZPtoT;

		APTemp[i][j][k].at     = (Dt - Ft/2.0); //max(max(-Ft , (Dt - Ft/2.0)), 0.0);
		APTemp[i][j][k+1].ab   = (Dt + Ft/2.0); //max(max( Ft , (Dt + Ft/2.0)), 0.0);
	
	}

	//for outlet boundary (top)

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = vec[i][j].size() - 2; k < vec[i][j].size() - 1; k++)
			{
				double Ft		= massFluxTop[i][j][k].value;
				double Dt		= (2*vec[i][j][k].visc*vec[i][j][k].St)/vec[i][j][k].DZPtoT;

				APTemp[i][j][k].at	= -Ft + Dt;
				
			}
		}
	
	}

	
	//for inlet boundary (bottom)
	
	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < 2; k++)
			{
				double Fb		=  massFluxTop[i][j][k-1].value;
				double Db		=  (2*vec[i][j][k].visc*vec[i][j][k].St)/vec[i][j][k].DZPtoT;

				APTemp[i][j][k].ab 	= Fb + Db;

			}
		}
	
	}


	//wall north

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = vec[i].size() - 2; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				double Fn 		= massFluxNorth[i][j][k].value;
				double Dn 		= (2*vec[i][j][k].visc*vec[i][j][k].Sn)/vec[i][j][k].DYPtoN;
	
				APTemp[i][j][k].an      =  max(max(-Fn , (Dn - Fn/2.0)), 0.0);

			}
		}

	}

	//Wall south

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<2; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				double Fs 		= massFluxNorth[i][j-1][k].value;
				double Ds 		= (2*vec[i][j][k].visc*vec[i][j][k].Sn)/vec[i][j][k].DYPtoN;
			
				APTemp[i][j][k].as   	=  max(max( Fs , (Ds + Fs/2.0)), 0.0);

			}
		}
	}
	//Wall east

	for(unsigned int i = vec.size() - 2; i<vec.size() - 1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				double Fe 		= massFluxEast[i][j][k].value;
				double De 		= (2*vec[i][j][k].visc*vec[i][j][k].Se)/vec[i][j][k].DXPtoE;

				APTemp[i][j][k].ae    	=  max(max(-Fe , (De - Fe/2.0)), 0.0);

			}
		}
	}

	//Wall west

	for(unsigned int i = 1; i < 2; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = 1; k < vec[i][j].size() - 1; k++)
			{
				double Fw 		= massFluxEast[i-1][j][k].value;
				double Dw 		= (2*vec[i][j][k].visc*vec[i][j][k].Se)/vec[i][j][k].DXPtoE;

				APTemp[i][j][k].aw     	=  max(max(Fw , (Dw + Fw/2.0)), 0.0);

			}
		}

	}
	return APTemp;

} // end convection- diffusive term


FiniteMatrix::finiteMat pressureGrad(Fields::vec3dField& vec, int& direction_) // receives the pressure field
{
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), vector<FiniteMatrix>(vec[0][0].size())));
	
	
	forAllInternal(vec)
	{
				
		//double DX = vec[i][j][k].X - vec[i-1][j][k].X;
		//double DY = vec[i][j][k].Y - vec[i][j-1][k].Y;
		//double DZ = vec[i][j][k].Z - vec[i][j][k-1].Z;

		

        	double pressureEastFace = (vec[i+1][j][k].value*vec[i][j][k].FXE) + (vec[i][j][k].value* (1.0-vec[i][j][k].FXE));
        	double pressureWestFace = (vec[i][j][k].value*vec[i-1][j][k].FXE) + (vec[i-1][j][k].value*(1.0-vec[i-1][j][k].FXE));


        	double pressureNorthFace = (vec[i][j+1][k].value*vec[i][j][k].FYN) + (vec[i][j][k].value*(1.0-vec[i][j][k].FYN));
        	double pressureSouthFace = (vec[i][j][k].value*vec[i][j-1][k].FYN) + (vec[i][j-1][k].value*(1.0-vec[i][j-1][k].FYN));

		double pressureTopFace    = (vec[i][j][k+1].value*vec[i][j][k].FZT) + (vec[i][j][k].value*(1.0-vec[i][j][k].FZT));
        	double pressureBottomFace = (vec[i][j][k].value*vec[i][j][k-1].FZT) + (vec[i][j][k-1].value*(1.0-vec[i][j][k-1].FZT));

        	double pressureEastGrad  = (pressureEastFace-pressureWestFace);
       		double pressureNorthGrad = (pressureNorthFace-pressureSouthFace);
		double pressureTopGrad   = (pressureTopFace - pressureBottomFace);

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
		
	return APTemp;
} // end pressureGrad



FiniteMatrix::finiteMat forceSource(Fields::vec3dField& vec, Solution& sol_, int& direction_) // receives the velocity fields to calculate the coriolis force and centrifugal force
{
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), vector<FiniteMatrix>(vec[0][0].size())));
	
	
	forAllInternal(vec)
	{
		if(direction_ == 1)  //for east direction (pass W velocity field )
        	{
            		APTemp[i][j][k].sf = -2*vec[i][j][k].density*vec[i][j][k].value*sol_.omega*vec[i][j][k].volume  
					     + vec[i][j][k].density*sol_.omega*sol_.omega*vec[i][j][k].Se*(vec[i][j][k].X*vec[i][j][k].X - vec[i-1][j][k].X*vec[i-1][j][k].X)/2;

       		}
		
        	else if(direction_ == 2) //for north direction
        	{
            		APTemp[i][j][k].sf = 0.0;
        	}

		else if(direction_ == 3) // for top direction (pass U velocity field)
		{
			APTemp[i][j][k].sf = 2*vec[i][j][k].density*vec[i][j][k].value*sol_.omega*vec[i][j][k].volume  
					     + vec[i][j][k].density*sol_.omega*sol_.omega*vec[i][j][k].St*(vec[i][j][k].Z*vec[i][j][k].Z - vec[i][j][k-1].Z*vec[i][j][k-1].Z)/2;
		}

		
        }
		
	return APTemp;
} // end forceSource


//  for time coeff 

FiniteMatrix::finiteMat timeCoeff(Fields::vec3dField& vec, Solution& sol_) // receives the velocity fields to calculate the coriolis force and centrifugal force
{
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), vector<FiniteMatrix>(vec[0][0].size())));
	
	
	forAllInternal(vec)
	{
		APTemp[i][j][k].apnot = sol_.density*vec[i][j][k].volume/sol_.dt;
		
        }
		
	return APTemp;
} // end timeCoeff



// for pressure correction **************its not massflux wall velocity****************
FiniteMatrix::finiteMat pressureCorrectionCoeff(Fields::vec3dField& vec, Fields::vec3dField& massFluxEast, Fields::vec3dField& massFluxNorth, Fields::vec3dField& massFluxTop)
{
	// recieves the massfluxes and velocity field ( only for distance, viscosity and surface area)
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), vector<FiniteMatrix>(vec[0][0].size())));

	forAllInternalUCVs(vec)
	{
		APTemp[i][j][k].ae = vec[i][j][k].density*vec[i][j][k].Se*massFluxEast[i][j][k].de;
		APTemp[i+1][j][k].aw = vec[i+1][j][k].density*vec[i+1][j][k].Se*massFluxEast[i][j][k].de;
	}

	forAllInternalVCVs(vec)
	{
		APTemp[i][j][k].an = vec[i][j][k].density*vec[i][j][k].Sn*massFluxNorth[i][j][k].dn;
		APTemp[i][j+1][k].as = vec[i][j+1][k].density*vec[i][j+1][k].Sn*massFluxNorth[i][j][k].dn;
	}

	forAllInternalWCVs(vec)
	{
		APTemp[i][j][k].at = vec[i][j][k].density*vec[i][j][k].St*massFluxTop[i][j][k].dt;
		APTemp[i][j][k+1].ab = vec[i][j][k+1].density*vec[i][j][k+1].St*massFluxTop[i][j][k].dt;
	}

	//for outlet boundary (top)

	for(unsigned int i = 1; i<vec.size()-1; i++)
	{
		for(unsigned int j = 1; j<vec[i].size() - 1; j++)
		{
			for(unsigned int  k = vec[i][j].size() - 2; k < vec[i][j].size() - 1; k++)
			{

				APTemp[i][j][k].at	= vec[i][j][k].density*vec[i][j][k].St*massFluxTop[i][j][k].dt;
				
			}
		}
	
	}

	
	//for inlet boundary (bottom)
	
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


	//wall north

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

	//Wall south

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
	//Wall east

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

	//Wall west

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

} // end pressure correction


FiniteMatrix::finiteMat pressureCorrectionSource(Fields::vec3dField& vec, Fields::vec3dField& massFluxEast, Fields::vec3dField& massFluxNorth, Fields::vec3dField& massFluxTop)
{
	FiniteMatrix::finiteMat APTemp(vec.size(), vector<vector<FiniteMatrix>>(vec[0].size(), vector<FiniteMatrix>(vec[0][0].size())));
	
	
	forAllInternal(vec)
	{
	
            	APTemp[i][j][k].sp = -massFluxEast[i][j][k].value + massFluxEast[i-1][j][k].value - massFluxNorth[i][j][k].value + massFluxNorth[i][j-1][k].value - massFluxTop[i][j][k].value + massFluxTop[i][j][k-1].value;

        }
		
	return APTemp;


}

} // end namespace fvm
