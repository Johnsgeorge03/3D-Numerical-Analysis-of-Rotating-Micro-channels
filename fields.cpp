#include "fields.hpp"


Fields::Fields()
:value(0.0)
{

}



Fields::Fields(int& NI_, int& NJ_, int& NK_)
:value(0.0), NI(NI_), NJ(NJ_), NK(NK_), NIM(NI-1), NJM(NJ-1), NKM(NK-1), de(0.0), dn(0.0), dt(0.0) //initializing some variables
{

}



Fields::~Fields()
{
}

//function to pass the grid and solution info the 3dvector of fields

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



void Fields::print3dmat(Fields::vec3dField& vec)
{
	for(unsigned int i = 0; i<vec.size(); i++)
	{
		for(unsigned int j = 0; j<vec[i].size(); j++)
		{
			for ( unsigned int k = 0; k<vec[i][j].size(); k++)
			{
				cout << vec[i][j][k].value << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
	}
}


void Fields::linearextrapolateCondition(Fields::vec3dField& vec, vector<double>& FXvec, vector<double>& FYvec, vector<double>& FZvec, string& wallname)
{
	if(wallname == "EAST")
	{
		forEastBoundary(vec)
		{
			vec[i][j][k].value = vec[i-1][j][k].value + (vec[i-1][j][k].value -
				 	   vec[i-2][j][k].value)*FXvec[i-1];
		}
	}

	else if(wallname == "WEST")
	{
		forWestBoundary(vec)
		{
			vec[i][j][k].value = vec[i+1][j][k].value + (vec[i+1][j][k].value -
				 	   vec[i+2][j][k].value)*FXvec[i+1];
		}
	}

	else if(wallname == "NORTH")
	{
		forNorthBoundary(vec)
		{

			vec[i][j][k].value = vec[i][j-1][k].value + (vec[i][j-1][k].value -
				 	   vec[i][j-2][k].value)*FYvec[j-1];

		}
	}

	else if(wallname == "SOUTH")
	{
		forSouthBoundary(vec)
		{

			vec[i][j][k].value = vec[i][j+1][k].value + (vec[i][j+1][k].value -
				 	   vec[i][j+2][k].value)*FYvec[j+1];
		}
	}

	else if(wallname == "TOP")
	{
		forTopBoundary(vec)
		{
			
			vec[i][j][k].value = vec[i][j][k-1].value + (vec[i][j][k-1].value -
				 	   vec[i][j][k-2].value)*FZvec[k-1];

		}

	}

	else if(wallname == "BOTTOM")
	{
		forBottomBoundary(vec)
		{
			
			vec[i][j][k].value = vec[i][j][k+1].value + (vec[i][j][k+1].value -
				 	   vec[i][j][k+2].value)*FZvec[k+1];
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
