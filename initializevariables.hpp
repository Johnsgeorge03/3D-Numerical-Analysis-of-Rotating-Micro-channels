// initialiazing geometry parameters

int NX = 30, NY = 30, NZ = 150;
double lengthX = 0.004, lengthY = 0.004, lengthZ = 0.7;
double beginX = 0.0, beginY = 0.0 , beginZ = 0.02;

Mesh mymesh_ (beginX, beginY, beginZ, NX, NY, NZ, lengthX, lengthY, lengthZ);


int NI , NJ, NK;
NI = mymesh_.getNI();
NJ = mymesh_.getNJ();
NK = mymesh_.getNK();

Solution sol;
Fields fieldOper(NI, NJ, NK);

// variables at the start of 
Fields::vec3dField U (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField V (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField W (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField P (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

// variables at the start of timestep or previous timestep
Fields::vec3dField UT (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField VT (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK))); 
Fields::vec3dField WT (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField PT (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

// cell centre velocity correction
Fields::vec3dField UC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField VC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField WC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

// pressure correction
Fields::vec3dField PC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

// wall velocity at the start of timestep or previous timestep
Fields::vec3dField UWT (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField VWT (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField WWT (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));

// face velocity and its correctios C stands for correction W for Wall
Fields::vec3dField UW (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField UWC (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField VW (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField VWC (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField WW (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));
Fields::vec3dField WWC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));

// massflux and its correction 
Fields::vec3dField massFluxE (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField massFluxN (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField massFluxT (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));
Fields::vec3dField massFluxEC (NI - 1, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
Fields::vec3dField massFluxNC (NI, Fields::vec2dField(NJ - 1, Fields::vec1dField(NK)));
Fields::vec3dField massFluxTC (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK - 1)));

// passing details about the mesh and solution parameter like density visc etc...
fieldOper.getGridInfoPassed(U, mymesh_, sol);
fieldOper.getGridInfoPassed(V, mymesh_, sol);
fieldOper.getGridInfoPassed(W, mymesh_, sol);

fieldOper.getGridInfoPassed(UT, mymesh_, sol); 
fieldOper.getGridInfoPassed(VT, mymesh_, sol);
fieldOper.getGridInfoPassed(WT, mymesh_, sol);
fieldOper.getGridInfoPassed(PT, mymesh_, sol);

fieldOper.getGridInfoPassed(P, mymesh_, sol);
fieldOper.getGridInfoPassed(PC, mymesh_, sol);
fieldOper.getGridInfoPassed(massFluxE, mymesh_, sol);
fieldOper.getGridInfoPassed(massFluxN, mymesh_, sol);
fieldOper.getGridInfoPassed(massFluxT, mymesh_, sol);

// used in setting boundary conditions
    string north  = "NORTH"; 
    string south  = "SOUTH";
    string east   = "EAST"; 
    string west   = "WEST";
    string top    = "TOP"; 
    string bottom = "BOTTOM";

// initial guess and inlet velocity
    double inletwvel  = 0.5;
    double u_guess    = 0.0;
    double v_guess    = 0.0;
    double w_guess    = 0.0;
    double massfluxwinlet = sol.density*inletwvel*U[1][1][1].St; // change the inlet velocity here
    double massfluxwguess = sol.density*w_guess*U[1][1][1].St;
    double massfluxuguess = sol.density*u_guess*U[1][1][1].Se;
    double massfluxvguess = sol.density*v_guess*U[1][1][1].Sn;

    //-------------------------------------BOUNDARY CONDITIONS-------------------------------------------//

    	//----------------- WALL ------------------//

	fieldOper.boundaryCondition(U, north, 0.0);
	fieldOper.boundaryCondition(U, south, 0.0);
	fieldOper.boundaryCondition(U, east, 0.0);
	fieldOper.boundaryCondition(U, west, 0.0);

	fieldOper.boundaryCondition(V, north, 0.0);
	fieldOper.boundaryCondition(V, south, 0.0);
	fieldOper.boundaryCondition(V, east, 0.0);
	fieldOper.boundaryCondition(V, west, 0.0);

	fieldOper.boundaryCondition(W, north, 0.0);
	fieldOper.boundaryCondition(W, south, 0.0);
	fieldOper.boundaryCondition(W, east, 0.0);
	fieldOper.boundaryCondition(W, west, 0.0);

	fieldOper.boundaryCondition(UT, north, 0.0);
	fieldOper.boundaryCondition(UT, south, 0.0);
	fieldOper.boundaryCondition(UT, east, 0.0);
	fieldOper.boundaryCondition(UT, west, 0.0);

	fieldOper.boundaryCondition(VT, north, 0.0);
	fieldOper.boundaryCondition(VT, south, 0.0);
	fieldOper.boundaryCondition(VT, east, 0.0);
	fieldOper.boundaryCondition(VT, west, 0.0);

	fieldOper.boundaryCondition(WT, north, 0.0);
	fieldOper.boundaryCondition(WT, south, 0.0);
	fieldOper.boundaryCondition(WT, east, 0.0);
	fieldOper.boundaryCondition(WT, west, 0.0);

	fieldOper.boundaryCondition(UW, north, 0.0);
	fieldOper.boundaryCondition(UW, south, 0.0);
	fieldOper.boundaryCondition(UW, east, 0.0);
	fieldOper.boundaryCondition(UW, west, 0.0);

	fieldOper.boundaryCondition(VW, north, 0.0);
	fieldOper.boundaryCondition(VW, south, 0.0);
	fieldOper.boundaryCondition(VW, east, 0.0);
	fieldOper.boundaryCondition(VW, west, 0.0);

	fieldOper.boundaryCondition(WW, north, 0.0);
	fieldOper.boundaryCondition(WW, south, 0.0);
	fieldOper.boundaryCondition(WW, east, 0.0);
	fieldOper.boundaryCondition(WW, west, 0.0);

	fieldOper.boundaryCondition(UWT, north, 0.0);
	fieldOper.boundaryCondition(UWT, south, 0.0);
	fieldOper.boundaryCondition(UWT, east, 0.0);
	fieldOper.boundaryCondition(UWT, west, 0.0);

	fieldOper.boundaryCondition(VWT, north, 0.0);
	fieldOper.boundaryCondition(VWT, south, 0.0);
	fieldOper.boundaryCondition(VWT, east, 0.0);
	fieldOper.boundaryCondition(VWT, west, 0.0);

	fieldOper.boundaryCondition(WWT, north, 0.0);
	fieldOper.boundaryCondition(WWT, south, 0.0);
	fieldOper.boundaryCondition(WWT, east, 0.0);
	fieldOper.boundaryCondition(WWT, west, 0.0);


	fieldOper.boundaryCondition(massFluxE, east, 0.0);
	fieldOper.boundaryCondition(massFluxE, west, 0.0);
	fieldOper.boundaryCondition(massFluxN, north, 0.0);
	fieldOper.boundaryCondition(massFluxN, south, 0.0);
	
	//-----------------END-----------------------//




	//------------BOTTOM OR INLET----------------//

	fieldOper.boundaryCondition(U, bottom, 0.0);
	fieldOper.boundaryCondition(V, bottom, 0.0);
	fieldOper.boundaryCondition(W, bottom, inletwvel);

	fieldOper.boundaryCondition(UT, bottom, 0.0);
	fieldOper.boundaryCondition(VT, bottom, 0.0);
	fieldOper.boundaryCondition(WT, bottom, inletwvel); 

	fieldOper.boundaryCondition(UW, bottom, 0.0);
	fieldOper.boundaryCondition(VW, bottom, 0.0);
	fieldOper.boundaryCondition(WW, bottom, inletwvel);

	fieldOper.boundaryCondition(UWT, bottom, 0.0);
	fieldOper.boundaryCondition(VWT, bottom, 0.0);
	fieldOper.boundaryCondition(WWT, bottom, inletwvel);

	fieldOper.boundaryCondition(massFluxT, bottom, massfluxwinlet);

	//------------------- END --------------------//





	//--------------- TOP OR OUTLET --------------//

	fieldOper.boundaryCondition(P, top, 0.0);
	fieldOper.boundaryCondition(PT, top, 0.0);

	//------------------- END --------------------//





	//----------INITIALIZE INTERNAL FIELD---------//
	
	fieldOper.initializeInternalField(U, u_guess);
	fieldOper.initializeInternalField(V, v_guess);
	fieldOper.initializeInternalField(W, w_guess);
	fieldOper.initializeInternalField(P, 0.0);

	// initialize initial time step zero as this is not guess

	fieldOper.initializeInternalField(UT, 0.0);
	fieldOper.initializeInternalField(VT, 0.0);
	fieldOper.initializeInternalField(WT, 0.0);
	fieldOper.initializeInternalField(PT, 0.0);

	fieldOper.initializeInternalField(UW, u_guess);
	fieldOper.initializeInternalField(VW, v_guess);
	fieldOper.initializeInternalField(WW, w_guess);

	fieldOper.initializeInternalField(UWT, 0.0);
	fieldOper.initializeInternalField(VWT, 0.0);
	fieldOper.initializeInternalField(WWT, 0.0);

	fieldOper.initializeInternalField(massFluxE, massfluxuguess);
	fieldOper.initializeInternalField(massFluxN, massfluxvguess);
	fieldOper.initializeInternalField(massFluxT, massfluxwguess);

	//------------------END-----------------------//

