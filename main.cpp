#include<iostream>
#include<cmath>
#include "mesh.hpp"
#include "solution.hpp"
#include "foralloperations.hpp"
#include "fields.hpp"
#include "finitematrix.hpp"
#include "equation.hpp"
#include "finitevolumeoperations.hpp"
#include "filewrite.hpp"
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
	#include "initializevariables.hpp"
    	#include "initializefinitematrixvar.hpp" 
	int timeiter = 0;
	int maxtimeit = 1; // max timestep

	//--------------------------------------- START OF TIME ITERATION -----------------------------------//

	for(timeiter = 0 ; timeiter< maxtimeit; timeiter++)
	{
	cout<<" -----------START OF TIMESTEP---------: " <<timeiter + 1<<"------------"<< endl;
	if(timeiter != 0)
	{
		U = UT;
		V = VT;
		W = WT;
		P = PT;
		UW = UWT;
		VW = VWT;
		WW = WWT;
	}
	double R2UMax = 0.0;
	double R2VMax = 0.0;
	double R2WMax = 0.0;
	double R2PMax = 0.0;
	double vel_tol = 1e-6;  //tolerance for solver - velocity
	double pres_tol = 1e-4; //tolerance for solver - pressure
	double tolerance = 1e-13; // tolerance for outer iteration
	int iter = 0;

	//--------------------------------------- START OF SIMPLE ITERATION ----------------------------------//

	for(iter = 0; iter< sol.maxit; iter++)
	{	

		cout<<" --------START OF OUTER ITERATION NO: "<<iter+1<<"--------------"<<endl;

		int eastdir  = 1;
		int northdir = 2; 
		int topdir   = 3;

		unsigned int Uiterations  = 10;
		unsigned int Viterations  = 10; 
		unsigned int Witerations  = 10;
		unsigned int pressureiter = 2000;
		

		
		//------------------------COEFFICIENT CALCULATION AND SOLVING MOMN EQN------------------------//

		Equation *UEqn = new Equation
		(
			fvm::convDiffusiveTerm(U, massFluxE, massFluxN, massFluxT) 
			+ fvm::pressureGrad(P, U, V, W, massFluxT, eastdir)//velocities are not used in the calc
			+ fvm::forceSource(W, sol, eastdir) // W is used in this calc
			+ fvm::timeCoeff(U, sol) // U is not used in the calc
 
		);
		UEqn->assembleEquation(massFluxE, massFluxN, massFluxT, U, V, W, iter); // computing AP
		UEqn->relax(U, UT); //under-relaxation
		cout<<" Solving U-momentum equation " <<endl;
		cout<<endl;
		Fields::vec3dField U_ (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
		U_ = U;

		if(iter == 0)
		{
			U_ = UEqn->sipSolver(U_, UT, sol, Uiterations, vel_tol);
		}
		else
		{
			U_ = UEqn->CGSTABvel(U_, UT, sol, Uiterations, vel_tol);
		}
		
		fieldOper.getGridInfoPassed(U_, mymesh_, sol);


		

		Equation *VEqn = new Equation
		(
			fvm::convDiffusiveTerm(V, massFluxE, massFluxN, massFluxT)
			+ fvm::pressureGrad(P, U, V, W, massFluxT, northdir)// velocities are not used in this calc
			+ fvm::forceSource(U, sol, northdir) // U is not used in calc
			+ fvm::timeCoeff(V, sol)  // V is not used in calc
		);
		VEqn->EqnName = "V-Eqn";
		VEqn->assembleEquation(massFluxE, massFluxN, massFluxT, U, V, W, iter);
		VEqn->relax(V, VT); // under-relaxation
		cout<<" Solving V-momentum equation "<<endl;
		cout<<endl;
		Fields::vec3dField V_ (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
		V_ = V;

		
		if(iter>0)
		{
			V_ = VEqn->CGSTABvel(V_, VT, sol, Viterations, vel_tol);
		}
		else
		{
			V_ = VEqn->sipSolver(V_, VT, sol, Viterations, vel_tol);
		}		
		fieldOper.getGridInfoPassed(V_, mymesh_, sol);




		Equation *WEqn = new Equation
		(
			fvm::convDiffusiveTerm(W, massFluxE, massFluxN, massFluxT)
			+ fvm::pressureGrad(P, U, V, W, massFluxT, topdir) // velocities are not used in the calc
			+ fvm::forceSource(U_, sol, topdir) // U_ is used here 
			+ fvm::timeCoeff(W, sol) // W is not used in calc
		);
		WEqn->EqnName = "W-Eqn";
		WEqn->assembleEquation(massFluxE, massFluxN, massFluxT, U, V, W, iter);
		WEqn->relax(W, WT);
		cout<< " Solving W-momentum equation "<< endl;
		cout<< endl;
		Fields::vec3dField W_ (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));
		W_ = W;


		if(iter > 0)
		{
			W_ = WEqn->CGSTABvel(W_, WT, sol, Witerations, vel_tol);
		}
		else
		{
			W_ = WEqn->sipSolver(W_, WT, sol,Witerations, vel_tol);
		}
		fieldOper.getGridInfoPassed(W_, mymesh_, sol);
		fieldOper.linearextrapolateCondition(W_, mymesh_.FX, mymesh_.FY, mymesh_.FZ, top);
		

		//-------------END OF COEFFICIENT CALCULATION AND SOLVING MOMN EQN ----------------------//

/*

		// TESTING
		fileWriter fileWriter3_;
		string filename3 = "zbeforecorrection";
		int time3 = iter;
		fileWriter3_.writeUVWP(filename3, time3, mymesh_,U_,V_,W_,P);
		// TESTING END
*/
	


		//--------------------- RHIE-CHOW MOMENTUM INTERPOLATION --------------------------------//

		UW = UEqn->momentumInterpolation(U_, P, U, UWT, UT, UW, VW, WW, sol, eastdir);
		VW = VEqn->momentumInterpolation(V_, P, V, VWT, VT, UW, VW, WW, sol, northdir);
		WW = WEqn->momentumInterpolation(W_, P, W, WWT, WT, UW, VW, WW, sol, topdir);
		fieldOper.boundaryCondition(WW, bottom, inletwvel);
		
		//------------------ END OF RHIE-CHOW MOMENTUM INTERPOLATION ----------------------------//	
		




		//------------------------- MASS FLUX CALCULATON ----------------------------------------//

		// additional vel field is sent to extract the density and face area//
		massFluxE = UEqn->massFluxCalculation(UW, U_, eastdir);
		massFluxN = VEqn->massFluxCalculation(VW, V_, northdir);
		massFluxT = WEqn->massFluxCalculation(WW, W_, topdir);
		fieldOper.boundaryCondition(massFluxT, bottom, massfluxwinlet);

		//----------------------- END OF MASS FLUX CALCULATION ----------------------------------//
	





		//-------------------------PRESSURE COEFFICIENT ASSEMBLY---------------------------------//

		cout<<" Calculating pressure coefficients " <<endl;
		cout<<endl;

		Equation *PEqn =  new Equation
		(
			fvm::pressureCorrectionCoeff(P, UW, VW, WW, UEqn->AP, VEqn->AP, WEqn->AP)
			+ fvm::pressureCorrectionSource(P, massFluxE, massFluxN, massFluxT)
		);

		PEqn->EqnName = "P-Eqn";
		PEqn->assemblePressureCorrectionEquation(WEqn->AP, P);
	//	PEqn->SOR = 1.0e-3;
		//-------------------- END OF PRESSURE COEFFICIENT ASSEMBLY -----------------------------//
		




		//--------------------INITIALIZING CORRECTION FIELD TO ZERO -----------------------------//

		fieldOper.initializeField(PC, 0.0);
		fieldOper.initializeField(UWC, 0.0);
		fieldOper.initializeField(VWC, 0.0);
		fieldOper.initializeField(WWC, 0.0);
		fieldOper.initializeField(UC, 0.0);
		fieldOper.initializeField(VC, 0.0);
		fieldOper.initializeField(WC, 0.0);
		fieldOper.initializeField(massFluxEC, 0.0);
		fieldOper.initializeField(massFluxNC, 0.0);
		fieldOper.initializeField(massFluxTC, 0.0);
		
		//------------------END OF INITIALIZING CORRECTION FIELD TO ZERO ------------------------//




	
		//--------------------------- SOLVING PRESSURE CORRECTION -------------------------------//
		
		cout<<" Solving pressure correction equation "<< endl;
		cout<<endl;
		PC = PEqn->CGSTAB(PC, sol, pressureiter, pres_tol);

/*
		if(iter > 0)
		{
			PC = PEqn->CGSTAB(PC, sol, pressureiter, pres_tol);
		}
		else
		{
			PC = PEqn->sipSolverPressure(PC, sol, pressureiter, pres_tol);
		}
*/


//		PC = PEqn->CGSTAB(PC, sol, pressureiter, pres_tol);
		fieldOper.getGridInfoPassed(PC, mymesh_, sol);
/*
		//----TESTING 
		cout<< " PC "<<" iteration "<<iter <<endl;
		fieldOper.print3dmat(PC);
		cout<<" end PC " <<endl;
		exit(0);
		//----END TESTING
*/

/*		
		PC = PEqn->ICCG(PC, sol, pressureiter, pres_tol);
		PC = PEqn->sipSolverPressure(PC, sol, pressureiter, pres_tol);
		PC = PEqn->solvePressure(PC, pressureiter);
*/		
		//---------------------------END OF SOLVING PRESSURE CORRECTION -------------------------//





		//--------------------------PRESSURE CORRECTION EXTRAPOLATION ---------------------------//

		fieldOper.linearextrapolateCondition(PC, mymesh_.FX, mymesh_.FY, mymesh_.FZ, east);
		fieldOper.linearextrapolateCondition(PC, mymesh_.FX, mymesh_.FY, mymesh_.FZ, west);
		fieldOper.linearextrapolateCondition(PC, mymesh_.FX, mymesh_.FY, mymesh_.FZ, north);
		fieldOper.linearextrapolateCondition(PC, mymesh_.FX, mymesh_.FY, mymesh_.FZ, south);
		fieldOper.linearextrapolateCondition(PC, mymesh_.FX, mymesh_.FY, mymesh_.FZ, bottom);

		//--------------------END OF PRESSURE CORRECTION EXTRAPOLATION --------------------------//





		//---------------------------------FIELD CORRECTION--------------------------------------//


		cout<< " Correcting Fields "<<endl;
		cout<<endl;


		UWC = UEqn->faceVelocityCorrection(PC, UW, VW, WW, eastdir);
		VWC = VEqn->faceVelocityCorrection(PC, UW, VW, WW, northdir);
		WWC = WEqn->faceVelocityCorrection(PC, UW, VW, WW, topdir);

		massFluxEC = UEqn->massFluxCorrection(PC, UWC, VWC, WWC, eastdir);
		massFluxNC = VEqn->massFluxCorrection(PC, UWC, VWC, WWC, northdir);
		massFluxTC = WEqn->massFluxCorrection(PC, UWC, VWC, WWC, topdir);


		UC = UEqn->cellVelocityCorrection(PC, eastdir);
		VC = VEqn->cellVelocityCorrection(PC, northdir); 
		WC = WEqn->cellVelocityCorrection(PC, topdir);

		//********* end of correction of velocity and massflux******************//


		//*********** adding the correction********************//
		forAllInternal(UW)
		{
			UW[i][j][k].value += UWC[i][j][k].value;
			massFluxE[i][j][k].value += massFluxEC[i][j][k].value;
		}	


		forAllInternal(VW)
		{
			VW[i][j][k].value += VWC[i][j][k].value;
			massFluxN[i][j][k].value += massFluxNC[i][j][k].value;
		}
		

		forAllInternal(WW)
		{
			WW[i][j][k].value += WWC[i][j][k].value;
			massFluxT[i][j][k].value += massFluxTC[i][j][k].value;
		}

		

		
/*	
		for( unsigned int i = 1; i<WW.size() - 1;i++)
		{
			for(unsigned int j = 1; j<WW[0].size() -1; j++)
			{
				for(unsigned int k = 1 ; k< WW[0][0].size(); k++)
				{
					WW[i][j][k].value += WWC[i][j][k].value;
					massFluxT[i][j][k].value += massFluxTC[i][j][k].value;
				}
			}
		}
*/
/*
		//----TESTING 
		cout<< " WW "<<" iteration "<<iter <<endl;
		fieldOper.print3dmat(WW);
		cout<<" end WW " <<endl;

		//----END TESTING
*/

		Fields::vec3dField P_ (NI, Fields::vec2dField(NJ, Fields::vec1dField(NK)));

		fieldOper.getGridInfoPassed(P_, mymesh_, sol);

		forAllInternal(P)
		{
			P_[i][j][k].value = P[i][j][k].value + sol.URFPressure*(PC[i][j][k].value);
			// P is the previous outer iteration pressure 
			U_[i][j][k].value = U_[i][j][k].value + UC[i][j][k].value;		
			// U_ in the rhs is obtained from solving momentum equation
			V_[i][j][k].value = V_[i][j][k].value + VC[i][j][k].value;
			W_[i][j][k].value = W_[i][j][k].value + WC[i][j][k].value;
		}


		
		//---------------------------------- END OF CORRECTION ------------------------------------//




		//--------------------------- PRESSURE & W EXTRAPOLATION ----------------------------------//

		fieldOper.linearextrapolateCondition(P_, mymesh_.FX, mymesh_.FY, mymesh_.FZ, east);
		fieldOper.linearextrapolateCondition(P_, mymesh_.FX, mymesh_.FY, mymesh_.FZ, west);
		fieldOper.linearextrapolateCondition(P_, mymesh_.FX, mymesh_.FY, mymesh_.FZ, north);
		fieldOper.linearextrapolateCondition(P_, mymesh_.FX, mymesh_.FY, mymesh_.FZ, south);
		fieldOper.linearextrapolateCondition(P_, mymesh_.FX, mymesh_.FY, mymesh_.FZ, bottom);

		fieldOper.linearextrapolateCondition(W_, mymesh_.FX, mymesh_.FY, mymesh_.FZ, top);

		//-------------------------END OF PRESSURE & W EXTRAPOLATION -------------------------------//


/*
		//---------------------CORRECTION OUTLET TO SATISFY CONTINUITY -----------------------------//

		
		double total_outlet_mass_flux = 0.0;
		double total_inlet_mass_flux = lengthX*lengthY*sol.density*inletwvel;

		forTopBoundary(W_)
		{
			total_outlet_mass_flux += W_[i][j][k].value;
		}

		double K = (total_inlet_mass_flux/total_outlet_mass_flux);
		
		cout<<endl;
		cout<< " Inlet / outlet " << K<< endl;
		if(total_outlet_mass_flux != 0)
		{
		forTopBoundary(W_)
		{
			W_[i][j][k].value = W_[i][j][k].value*K;
			W_[i][j][k-1].value = W_[i][j][k-1].value*K;

		}
		}
		//---------------------------END OF CONTINUITY CORRECTION AT OUTLET------------------------//
*/

		

		//---------------------- COPYING OUTLET CELL CENTER VELOCITY TO FACE -----------------------//

		forTopBoundary(WW)
		{
			WW[i][j][k].value = W_[i][j][k+1].value;
			massFluxT[i][j][k].value = sol.density*W_[i][j][k].St*W_[i][j][k+1].value;
		}


		//----------------- END OF COPYING OUTLET CELL CENTER VELOCITY TO FACE ---------------------//







		//------------------WRITING FILES AND UPDATING VALUES FOR NEXT ITERATION ------------------//

		P = P_;
		U = U_;
		V = V_;
		W = W_;
/*
		forAllInternal(U)
		{
			U[i][j][k].value = U[i][j][k].value + UEqn->URF*(U_[i][j][k].value - U[i][j][k].value);
			V[i][j][k].value = V[i][j][k].value + VEqn->URF*(V_[i][j][k].value - V[i][j][k].value);
			W[i][j][k].value = W[i][j][k].value + WEqn->URF*(W_[i][j][k].value - W[i][j][k].value);
		}
		fieldOper.linearextrapolateCondition(W, mymesh_.FX, mymesh_.FY, mymesh_.FZ, top);
*/

		fileWriter fileWriter2_;
		string filename = "zsimpleiteration";
		int time1 = iter;
		fileWriter2_.writeUVWP(filename, time1, mymesh_,U,V,W,P);
		
		//-----------------END OF WRITING FILES AND UPDATING FILES FOR NEXT ITERATION--------------//






		//--------------------------------- TEST FOR CONVERGENCE-----------------------------------//

		if(UEqn->RESOR >= R2UMax)
		{
			R2UMax = UEqn->RESOR;
		}

		if(VEqn->RESOR >= R2VMax)
		{
			R2VMax = VEqn->RESOR;
		}

		if(WEqn->RESOR >= R2WMax)
		{
			R2WMax = WEqn->RESOR;
		}

		if(PEqn->RESOR >= R2PMax)
		{
			R2PMax = PEqn->RESOR;
		}
		
	   	cout<< " IN OUTER ITERATION "<<iter+1<<" R2UNorm : " <<UEqn->RESOR/R2UMax<<" R2VNorm : "<< 
		VEqn->RESOR/R2VMax << " R2WNorm : "<< WEqn->RESOR/R2WMax<<" R2PNorm : "<<
		PEqn->RESOR/R2PMax<<endl;
	    	cout<<endl;

		if(PEqn->RESOR <= tolerance)
		{
			fileWriter fileWriter_;
			string filename = "zpipeflow_rotation";
			int time1 = timeiter;
			fileWriter_.writeUVWP(filename, time1, mymesh_, U, V, W, P);

			cout<< "-------*****----------RESULTS CONVERGED FOR OUTER-ITERATION "<< iter + 1 <<
				" FOR TIMESTEP " << timeiter + 1 <<"-----------**********-------------"<<endl;
			UT = U;
			VT = V;
			WT = W;
			PT = P;
			UWT = UW;
			VWT = VW;
			WWT = WW;
			delete UEqn;
			delete VEqn;
			delete WEqn;
			delete PEqn;
			break;
		}


		//------------------------------- END OF TEST FOR CONVERGENCE-------------------------------//



		delete UEqn;
		delete VEqn;
		delete WEqn;
		delete PEqn;

		cout<<" -------------END OF OUTER ITERATION NO: "<<iter+1<<"------------"<<endl;
		cout<<endl;
		cout<<endl;


	} 
	//------------------------------------ END OF SIMPLE ITERATION --------------------------------------//
	cout<<"------------END OF TIMESTEP: " << timeiter + 1<<endl;
	cout<<endl;
	cout<<endl;
	cout<<endl;
	} 
	//----------------------------------    END OF TIME ITERATION ---------------------------------------//
				
	
    return 0;

}
/*		cout<<" MassfluxE "<<endl;
		fieldOper.print3dmat(massFluxE);
		cout<<" MassfluxN "<<endl;
		fieldOper.print3dmat(massFluxN);
		cout<<" MassfluxT "<<endl;
		fieldOper.print3dmat(massFluxT);
*/
/*		cout<<" AP "<<endl;
		finiteObj.print3dmat(PEqn->AP);
		cout<<" AE "<<endl;
		finiteObj.print3dmat(PEqn->AE);
		cout<<" AW "<<endl;
		finiteObj.print3dmat(PEqn->AW);
		cout<<" AN "<<endl;
		finiteObj.print3dmat(PEqn->AN);
		cout<<" AS "<<endl;
		finiteObj.print3dmat(PEqn->AS);
		cout<<" AT "<<endl;
		finiteObj.print3dmat(PEqn->AT);
		cout<<" AB "<<endl;
		finiteObj.print3dmat(PEqn->AB);
		cout<<" SP "<<endl;
		finiteObj.print3dmat(PEqn->SP);
		exit(0);
*/

