#include "filewrite.hpp"
#include <iostream>
#include <string>

fileWriter::fileWriter()
{

}

fileWriter::~fileWriter()
{

}


void fileWriter::writeUVWP(string& name, int time_, Mesh& Mesh_, 
	Fields::vec3dField& Utemp, Fields::vec3dField& Vtemp, Fields::vec3dField& Wtemp,
	Fields::vec3dField& Ptemp)
{
	string myfileType	= ".dat";
	int ttime		= time_+1;
	string timename;
	ostringstream temp;
	temp<<ttime;
	timename 		= temp.str();
	string newnames_ 	= name;
	string newfilename 	= newnames_.append(timename);
	string name2 		= newfilename.append(myfileType);
	char cstr[name2.size() + 2];
	strcpy(cstr, newfilename.c_str());
	//   cout << "file name " << cstr << endl;
 	FILE *outfile;
	outfile 		= fopen( cstr,"w+t");

        int NXtemp 		= Utemp.size();
        int NYtemp 		= Utemp[0].size();
	int NZtemp		= Utemp[0][0].size();

	fprintf(outfile, "VARIABLES=\t\"X\"\t, \t\"Y\"\t, \t\"Z\"\t, \t\"U\"\t, \t\"V\"\t, \t\"W\"\t,\t\"P\"\n");
	fprintf(outfile, "ZONE  F=POINT\n");
	fprintf(outfile, "I=%d, J=%d, K=%d\n", NXtemp, NYtemp, NZtemp);
    	double xpos, ypos, zpos, UU, VV, WW, PP;


    	for ( unsigned int i = 0 ; i <Utemp.size() ; i++ )
    	{
    		for ( unsigned int j = 0 ;j<Utemp[0].size() ;j++ )
    		{
			for( unsigned int k = 0; k<Utemp[0][0].size(); k++)
			{
				xpos = Mesh_.XC[i];
				ypos = Mesh_.YC[j];
				zpos = Mesh_.ZC[k];
       				UU = Utemp[i][j][k].value;
        			VV = Vtemp[i][j][k].value;
				WW = Wtemp[i][j][k].value;
       				PP = Ptemp[i][j][k].value;
            //%5.8lf\t
			fprintf(outfile, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n",
				xpos, ypos, zpos, UU, VV, WW, PP);
    			}
		}

	}

	fclose(outfile);

}

