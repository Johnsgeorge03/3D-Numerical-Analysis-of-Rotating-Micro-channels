int NX = 20, NY = 20, NZ = 60;
double lengthX = 2e-5, lengthY = 3e-5, lengthZ = 9e-5;
double beginX = 0.0, beginY = 0.0 , beginZ = 0;

Mesh mymesh_ (beginX, beginY, beginZ, NX, NY, NZ, lengthX, lengthY, lengthZ);


int NI , NJ, NK;
NI = mymesh_.getNI();
NJ = mymesh_.getNJ();
NK = mymesh_.getNK();

Solution sol;
Fields fieldOper(NI, NJ, NK);

Fields::vec3dField U (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField V (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField W (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField P (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField UW (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField UWC (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField VW (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField VWC (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField WW (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));
Fields::vec3dField WWC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));
Fields::vec3dField PC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));  //pressure correction
Fields::vec3dField DPX (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField DPY (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK))); //pressure gradient in y direction
Fields::vec3dField DPZ (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField massFluxE (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField massFluxN (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField massFluxT (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));
Fields::vec3dField massFluxEC (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField massFluxNC (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField massFluxTC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));


fieldOper.getGridInfoPassed(U, mymesh_, sol);
fieldOper.getGridInfoPassed(V, mymesh_, sol);
fieldOper.getGridInfoPassed(W, mymesh_, sol);
fieldOper.getGridInfoPassed(P, mymesh_, sol);
fieldOper.getGridInfoPassed(PC, mymesh_, sol);
fieldOper.getGridInfoPassed(DPX, mymesh_, sol);
fieldOper.getGridInfoPassed(DPY, mymesh_, sol);
fieldOper.getGridInfoPassed(DPZ, mymesh_, sol);
fieldOper.getGridInfoPassed(massFluxE, mymesh_, sol);
fieldOper.getGridInfoPassed(massFluxN, mymesh_, sol);
fieldOper.getGridInfoPassed(massFluxT, mymesh_, sol);




    cout<<"***************************"<< endl;
   // cout<<"X"<<endl;
   // for( unsigned int i =0; i<mymesh_.X.size(); i++)
   // {
   //     cout<< mymesh_.X[i]<<" ";
   // }
   // cout<<"Y"<<endl;
   // 
   // for(unsigned int i = 0; i< mymesh_.Y.size(); i++)
   // {
   //     cout<< mymesh_.Y[i]<<" ";
   // }
   // cout<<"Z"<<endl;

   // for( unsigned int i =0; i<mymesh_.Z.size(); i++)
   // {
   //     cout<< mymesh_.Z[i]<<" ";
   // }
   // cout<<"XC"<<endl;

   // for( unsigned int i =0; i<mymesh_.XC.size(); i++)
   // {
   //     cout<< mymesh_.XC[i]<<" ";
   // }
   // cout<<"YC"<<endl;

   // for( unsigned int i =0; i<mymesh_.YC.size(); i++)
   // {
   //     cout<< mymesh_.YC[i]<<" ";
   // }
   // cout<<"ZC"<<endl;
   // for( unsigned int i =0; i<mymesh_.ZC.size(); i++)
   // {
   //     cout<< mymesh_.ZC[i]<<" ";
   // }
   // cout<<"FX"<<endl;
   // for( unsigned int i =0; i<mymesh_.FX.size(); i++)
   // {
   //     cout<< mymesh_.FX[i]<<" ";
   // }
   // cout<<"FY"<<endl;
   // for( unsigned int i =0; i<mymesh_.FY.size(); i++)
   // {
   //     cout<< mymesh_.FY[i]<<" ";
   // }
   // cout<<"FZ"<<endl;
   // for( unsigned int i =0; i<mymesh_.FZ.size(); i++)
   // {
   //     cout<< mymesh_.FZ[i]<<" ";
   // }
   // cout<<"U"<<endl;
    string north = "NORTH"; string south = "SOUTH";
    string east = "EAST"; string west = "WEST";
    string top = "TOP"; string bottom = "BOTTOM";
   // fieldOper.boundaryCondition(U,north, 1.0);
    //fieldOper.boundaryCondition(U,south, 2.0);
   // fieldOper.boundaryCondition(U,east, 3.0);
   // fieldOper.boundaryCondition(U,west, 4.0);
   // fieldOper.boundaryCondition(U,top, 5.0);
   // fieldOper.boundaryCondition(U,bottom, 6.0);
   // fieldOper.print3dmat(U);
    //fieldOper.initializeInternalField(P, 0.2);
   // fieldOper.copyInternalField(P, PP);
   // fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, north);
   // fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, south);
    //fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, east);
    //fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, west);
    //fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, top);
    //fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, bottom);

    //FiniteMatrix::finiteMat AE2(AE);
    //forAll(AE)
    //{
	//AE[i][j][k].value = 2.0;
	//AE2[i][j][k].value = 3.0;
    //}
    //const double two = 2.0;
    //FiniteMatrix::finiteMat AE3(AE2&&AE);
    
    //finiteObj.print3dmat(AE3);
                              
    //fieldOper.print3dmat(P);
    

    // Equations
	fieldOper.boundaryCondition(U, north, 0.0);
	fieldOper.boundaryCondition(U, south, 0.0);
	fieldOper.boundaryCondition(U, east, 0.0);
	fieldOper.boundaryCondition(U, west, 0.0);
	fieldOper.boundaryCondition(U, bottom, 0.0);

	fieldOper.boundaryCondition(V, north, 0.0);
	fieldOper.boundaryCondition(V, south, 0.0);
	fieldOper.boundaryCondition(V, east, 0.0);
	fieldOper.boundaryCondition(V, west, 0.0);
	fieldOper.boundaryCondition(V, bottom, 0.0);

	fieldOper.boundaryCondition(W, north, 0.0);
	fieldOper.boundaryCondition(W, south, 0.0);
	fieldOper.boundaryCondition(W, east, 0.0);
	fieldOper.boundaryCondition(W, west, 0.0);
	fieldOper.boundaryCondition(W, bottom, 0.1);


	fieldOper.boundaryCondition(UW, north, 0.0);
	fieldOper.boundaryCondition(UW, south, 0.0);
	fieldOper.boundaryCondition(UW, east, 0.0);
	fieldOper.boundaryCondition(UW, west, 0.0);
	fieldOper.boundaryCondition(UW, bottom, 0.0);

	fieldOper.boundaryCondition(VW, north, 0.0);
	fieldOper.boundaryCondition(VW, south, 0.0);
	fieldOper.boundaryCondition(VW, east, 0.0);
	fieldOper.boundaryCondition(VW, west, 0.0);
	fieldOper.boundaryCondition(VW, bottom, 0.0);

	fieldOper.boundaryCondition(WW, north, 0.0);
	fieldOper.boundaryCondition(WW, south, 0.0);
	fieldOper.boundaryCondition(WW, east, 0.0);
	fieldOper.boundaryCondition(WW, west, 0.0);
	fieldOper.boundaryCondition(WW, bottom, 0.1);

	fieldOper.boundaryCondition(UWC, north, 0.0);
	fieldOper.boundaryCondition(UWC, south, 0.0);
	fieldOper.boundaryCondition(UWC, east, 0.0);
	fieldOper.boundaryCondition(UWC, west, 0.0);
	fieldOper.boundaryCondition(UWC, bottom, 0.0);

	fieldOper.boundaryCondition(VWC, north, 0.0);
	fieldOper.boundaryCondition(VWC, south, 0.0);
	fieldOper.boundaryCondition(VWC, east, 0.0);
	fieldOper.boundaryCondition(VWC, west, 0.0);
	fieldOper.boundaryCondition(VWC, bottom, 0.0);

	fieldOper.boundaryCondition(WWC, north, 0.0);
	fieldOper.boundaryCondition(WWC, south, 0.0);
	fieldOper.boundaryCondition(WWC, east, 0.0);
	fieldOper.boundaryCondition(WWC, west, 0.0);
	fieldOper.boundaryCondition(WWC, bottom, 0.0);

	fieldOper.boundaryCondition(PC, top, 0.0);
	fieldOper.boundaryCondition(P, top, 0.0);

	fieldOper.boundaryCondition(massFluxE, east, 0.0);
	fieldOper.boundaryCondition(massFluxE, west, 0.0);
	fieldOper.boundaryCondition(massFluxN, north, 0.0);
	fieldOper.boundaryCondition(massFluxN, south, 0.0);
	//fieldOper.print3dmat(U);














