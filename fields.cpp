#include "fields.hpp"


Fields::Fields()
:value(0.0)
{

}



Fields::Fields(int& NI_, int& NJ_, int& NK_)
:value(0.0), NI(NI_), NJ(NJ_), NK(NK_), NIM(NI-1), NJM(NJ-1), NKM(NK-1), de(0.0), dn(0.0), dt(0.0) 
//initializing some variables
{

}



Fields::~Fields()
{
}


//function to pass the grid and solution info the 3dvector of fields
void Fields::copyAllField(Fields::vec3dField& from, Fields::vec3dField& to)
{

	forAll(from)
	{
		to[i][j][k].X   	= from[i][j][k].X;
		to[i][j][k].Y   	= from[i][j][k].Y;
		to[i][j][k].Z   	= from[i][j][k].Z;
      	 	to[i][j][k].XC  	= from[i][j][k].XC;
        	to[i][j][k].YC 		= from[i][j][k].YC;
        	to[i][j][k].ZC  	= from[i][j][k].ZC;
        	to[i][j][k].FXE 	= from[i][j][k].FXE;
        	to[i][j][k].FYN 	= from[i][j][k].FYN;
        	to[i][j][k].FZT 	= from[i][j][k].FZT;
        	to[i][j][k].FXP 	= from[i][j][k].FXP;
        	to[i][j][k].FYP 	= from[i][j][k].FYP;
        	to[i][j][k].FZP 	= from[i][j][k].FXP;
        	to[i][j][k].visc 	= from[i][j][k].visc;
       		to[i][j][k].density 	= from[i][j][k].density;
		to[i][j][k].value	= from[i][j][k].value;
		to[i][j][k].de		= from[i][j][k].de;
		to[i][j][k].dn		= from[i][j][k].dn;
		to[i][j][k].dt		= from[i][j][k].dt;
	}

	forAllInternal(from)
	{
		to[i][j][k].DXPtoE 	= from[i][j][k].DXPtoE;
		to[i][j][k].DYPtoN   	= from[i][j][k].DYPtoN;
		to[i][j][k].DZPtoT	= from[i][j][k].DZPtoT;
      	 	to[i][j][k].Se  	= from[i][j][k].Se;
        	to[i][j][k].Sn 		= from[i][j][k].Sn;
        	to[i][j][k].St  	= from[i][j][k].St;
        	to[i][j][k].volume 	= from[i][j][k].volume;
   
	}
}





void Fields::getGridInfoPassed(Fields::vec3dField& f, Mesh& mesh_, Solution& sol_)
{
    forAll(f)
    {
       f[i][j][k].X   	= mesh_.X[i];
       f[i][j][k].Y   	= mesh_.Y[j];
       f[i][j][k].Z   	= mesh_.Z[k];
       f[i][j][k].XC  	= mesh_.XC[i];
       f[i][j][k].YC 	= mesh_.YC[j];
       f[i][j][k].ZC  	= mesh_.ZC[k];
       f[i][j][k].FXE 	= mesh_.FX[i];
       f[i][j][k].FYN 	= mesh_.FY[j];
       f[i][j][k].FZT 	= mesh_.FZ[k];
       f[i][j][k].FXP 	= 1.0 - mesh_.FX[i];
       f[i][j][k].FYP 	= 1.0 - mesh_.FY[j];
       f[i][j][k].FZP 	= 1.0 - mesh_.FZ[k];
       f[i][j][k].visc 	= sol_.visc;
       f[i][j][k].density = sol_.density;
    }

    forAllInternal(f)
    {
        f[i][j][k].DXPtoE = mesh_.XC[i+1] - mesh_.XC[i];
        f[i][j][k].DYPtoN = mesh_.YC[j+1] - mesh_.YC[j];
        f[i][j][k].DZPtoT = mesh_.ZC[k+1] - mesh_.ZC[k];
        f[i][j][k].Se     = (mesh_.Y[j] - mesh_.Y[j-1])*(mesh_.Z[k] - mesh_.Z[k-1]);
        f[i][j][k].Sn     = (mesh_.X[i] - mesh_.X[i-1])*(mesh_.Z[k] - mesh_.Z[k-1]);    
        f[i][j][k].St     = (mesh_.Y[j] - mesh_.Y[j-1])*(mesh_.X[i] - mesh_.X[i-1]);
        f[i][j][k].volume = (mesh_.Y[j] - mesh_.Y[j-1])*(mesh_.Z[k] - mesh_.Z[k-1])*(mesh_.X[i] - mesh_.X[i-1]);

    }

}





void Fields::setVectorFieldGridFeatures()
{
}




//initialize the entire field

void Fields::initializeField(Fields::vec3dField& vec, double val)
{
	forAll(vec)
	{
		vec[i][j][k].value = val;
	}
}





//initialize the internal field only

void Fields::initializeInternalField(Fields::vec3dField& vec, double val)
{
	forAllInternal(vec)
	{
		vec[i][j][k].value = val;
	}
}





//to copy the internal field only
void Fields::copyInternalField(Fields::vec3dField& from, Fields::vec3dField& to)
{
	forAllInternal(from)
	{
		to[i][j][k].value = from[i][j][k].value;
	}
}





void Fields::boundaryCondition(Fields::vec3dField& vec,string& wallname, double bvalue)
{
	if(wallname == "EAST")
	{
		forEastBoundary(vec)
		{
			vec[i][j][k].value = bvalue;
		}
	}

	else if(wallname == "WEST")
	{
		forWestBoundary(vec)
		{
			vec[i][j][k].value = bvalue;
		}
	}

	else if(wallname == "NORTH")
	{
		forNorthBoundary(vec)
		{
			vec[i][j][k].value = bvalue;
		}
	}

	else if(wallname == "SOUTH")
	{
		forSouthBoundary(vec)
		{
			vec[i][j][k].value = bvalue;
		}
	}

	else if(wallname == "TOP")
	{
		forTopBoundary(vec)
		{
			vec[i][j][k].value = bvalue;
		}

	}

	else if(wallname == "BOTTOM")
	{
		forBottomBoundary(vec)
		{
			vec[i][j][k].value = bvalue;
		}
	}
}





void Fields::print3dmat(const Fields::vec3dField& vec)
{
	for(unsigned int i = 0; i<vec.size(); i++)
	{
		for(unsigned int j = 0; j<vec[i].size(); j++)
		{
			for ( unsigned int k = 0; k<vec[i][j].size(); k++)
			{
				cout << vec[i][j][k].value<< " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
	}
}





void Fields::linearextrapolateCondition(Fields::vec3dField& vec, vector<double>& FXvec, 
			vector<double>& FYvec, vector<double>& FZvec, string& wallname)
{
	if(wallname == "EAST")
	{
		forEastBoundary(vec)
		{
			vec[i][j][k].value = (3*vec[i-1][j][k].value - vec[i-2][j][k].value)*0.5;
		}
	}

	else if(wallname == "WEST")
	{
		forWestBoundary(vec)
		{
			vec[i][j][k].value = (3*vec[i+1][j][k].value - vec[i+2][j][k].value)*0.5;
		}
	}

	else if(wallname == "NORTH")
	{
		forNorthBoundary(vec)
		{

			vec[i][j][k].value = (3*vec[i][j-1][k].value - vec[i][j-2][k].value)*0.5;

		}
	}

	else if(wallname == "SOUTH")
	{
		forSouthBoundary(vec)
		{

			vec[i][j][k].value = (3*vec[i][j+1][k].value - vec[i][j+2][k].value)*0.5;
		}
	}

	else if(wallname == "TOP")
	{
		forTopBoundary(vec)
		{
			
			vec[i][j][k].value = vec[i][j][k-1].value;

		}

	}

	else if(wallname == "BOTTOM")
	{
		forBottomBoundary(vec)
		{
			
			vec[i][j][k].value = (3*vec[i][j][k+1].value - vec[i][j][k+2].value)*0.5;
		}
	}
}




void Fields::copyOutletVelocity(Fields::vec3dField& vec)
{
	forTopBoundary(vec)
	{
		vec[i][j][k].value = vec[i][j][k-1].value;
	}

}


