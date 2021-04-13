#include<iostream>
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
	
	int iter = 0;
	int maxit = 100; // timestep iteration
	for(iter = 0; iter< maxit; iter++)
	{	

		cout<<"start of iteration "<<iter<<endl;

		int eastdir = 1;
		int northdir = 2; 
		int topdir = 3;
		unsigned int iterations = 10;  // gauss - siedel iteration 
		unsigned int pressureiter = 30;

		
		
		//************Coefficient calculation**********************//
		Equation *UEqn = new Equation
		(
			fvm::convDiffusiveTerm(U, massFluxE, massFluxN, massFluxT) 
			+ fvm::pressureGrad(P, U, V, W, massFluxT, eastdir)
			+ fvm::forceSource(W, sol, eastdir)
			+ fvm::timeCoeff(U, sol)

		);
		
		UEqn->assembleEquation(massFluxE, massFluxN, massFluxT, U, V, W, iter); // computing AP
		
		Equation *VEqn = new Equation
		(
			fvm::convDiffusiveTerm(V, massFluxE, massFluxN, massFluxT)
			+ fvm::pressureGrad(P, U, V, W, massFluxT, northdir)
			+ fvm::forceSource(U, sol, northdir) // U is not used in calc
			+ fvm::timeCoeff(V, sol)
		);

		VEqn->assembleEquation(massFluxE, massFluxN, massFluxT, U, V, W, iter);

		Equation *WEqn = new Equation
		(
			fvm::convDiffusiveTerm(W, massFluxE, massFluxN, massFluxT)
			+ fvm::pressureGrad(P, U, V, W, massFluxT, topdir)
			+ fvm::forceSource(U, sol, topdir)
			+ fvm::timeCoeff(W, sol)
		);

		WEqn->assembleEquation(massFluxE, massFluxN, massFluxT, U, V, W, iter);
		//**************end of coefficient calculation********************//
//		cout<< "Wonly"<<endl;
//		fieldOper.print3dmat(W);
//		cout<<" Wonly end " <<endl;	
//		
//		cout<<" AP " <<endl;
//		finiteObj.print3dmat(UEqn->AP);
//		cout<<" AP end "<<endl;
//		
//		cout<<" AE " <<endl;
//		finiteObj.print3dmat(UEqn->AE);
//		cout<<" AE end "<<endl;
//		
//		cout<<" AW " <<endl;
//		finiteObj.print3dmat(UEqn->AW);
//		cout<<" AW end "<<endl;
//		
//		cout<<" AN " <<endl;
//		finiteObj.print3dmat(UEqn->AN);
//		cout<<" AN end "<<endl;
//		
//		cout<<" AS " <<endl;
//		finiteObj.print3dmat(UEqn->AS);
//		cout<<" AS end "<<endl;
//		
//		cout<<" AT " <<endl;
//		finiteObj.print3dmat(UEqn->AT);
//		cout<<" AT end "<<endl;
//		
//		cout<<" SP " <<endl;
//		finiteObj.print3dmat(UEqn->SP);
//		cout<<" SP end "<<endl;
//
//		cout<<" SF " <<endl;
//		finiteObj.print3dmat(WEqn->SF);
//		cout<<" SF end "<<endl;
//
//		cout<<" APNot " <<endl;
//		finiteObj.print3dmat(UEqn->APNot);
//		cout<<" APNot end "<<endl;
//
//		cout<<"rAP " <<endl;
//		finiteObj.print3dmat(UEqn->rAP);
//		cout<<"rAP end"<<endl;
		
//		cout<<"difference "<<endl;
//		finiteObj.print3dmat(UEqn->APU);
//		cout<<" end difference" <<endl;

		//**************Solving the equation that satisfy only momemtum equation**************//
		U_ = UEqn->solveVelocity(U, iterations);  // Gauss-Siedel iterations
		V_ = VEqn->solveVelocity(V, iterations);
		W_ = WEqn->solveVelocity(W, iterations);
		fieldOper.copyOutletVelocity(U_);
		fieldOper.copyOutletVelocity(V_);
		fieldOper.copyOutletVelocity(W_);

		fieldOper.getGridInfoPassed(U_, mymesh_, sol);
		fieldOper.getGridInfoPassed(V_, mymesh_, sol);
		fieldOper.getGridInfoPassed(W_, mymesh_, sol);

		fieldOper.boundaryCondition(U_, north, 0.0);
		fieldOper.boundaryCondition(U_, south, 0.0);
		fieldOper.boundaryCondition(U_, east, 0.0);
		fieldOper.boundaryCondition(U_, west, 0.0);
		fieldOper.boundaryCondition(U_, bottom, 0.0);

		fieldOper.boundaryCondition(V_, north, 0.0);
		fieldOper.boundaryCondition(V_, south, 0.0);
		fieldOper.boundaryCondition(V_, east, 0.0);
		fieldOper.boundaryCondition(V_, west, 0.0);
		fieldOper.boundaryCondition(V_, bottom, 0.0);

		fieldOper.boundaryCondition(W_, north, 0.0);
		fieldOper.boundaryCondition(W_, south, 0.0);
		fieldOper.boundaryCondition(W_, east, 0.0);
		fieldOper.boundaryCondition(W_, west, 0.0);
		fieldOper.boundaryCondition(W_, bottom, 0.01);

		//***********************end of solving partial equation******************************//
		
//		cout<< "U_"<<endl;
//		fieldOper.print3dmat(U_);
//		cout<<" U_ end" <<endl;	



		//****************momentum interpolation****************************//
		UW = UEqn->momentumInterpolation(U_, P, eastdir);
		VW = VEqn->momentumInterpolation(V_, P, northdir);
		WW = WEqn->momentumInterpolation(W_, P, topdir);
		
		//***************end of momentum interpolation**********************//	
		
//		cout<< "UW"<<endl;
//		fieldOper.print3dmat(UW);
//		cout<<" UW end" <<endl;	



		//****************mass flux calculation*****************************//
		// additional vel field is sent to extract the density and face area//
		massFluxE = UEqn->massFluxCalculation(UW, U, eastdir);
		massFluxN = VEqn->massFluxCalculation(VW, V, northdir);
		massFluxT = WEqn->massFluxCalculation(WW, W, topdir);
		//**************end of mass flux calculation***********************//
//		cout<< "massFluxE"<<endl;
//		fieldOper.print3dmat(massFluxE);
//		cout<<" massFluxE end" <<endl;	
	



		//***************pressure coefficient assembly*********************//
		Equation *PEqn =  new Equation
		(
			fvm::pressureCorrectionCoeff(P, UW, VW, WW)
			+ fvm::pressureCorrectionSource(P, massFluxE, massFluxN, massFluxT)
		);

		PEqn->assemblePressureEquation();
		//**************end of pressure coefficient assembly*******************//





		//*****************pressure correction solve************************//
		PC =  PEqn->solvePressure(PC, pressureiter);
		//*************end of pressure correction solve*********************//

//		cout<<" PC "<<endl;
//		fieldOper.print3dmat(PC);
//		cout<<" PC end "<<endl;	


			
		




		//*******************correction of velocity and massflux************//
		UWC = UEqn->faceVelocityCorrection(PC, UW, VW, WW, eastdir);
		VWC = VEqn->faceVelocityCorrection(PC, UW, VW, WW, northdir);
		WWC = WEqn->faceVelocityCorrection(PC, UW, VW, WW, topdir);

		massFluxEC = UEqn->massFluxCorrection(PC, UWC, VWC, WWC, eastdir);
		massFluxNC = VEqn->massFluxCorrection(PC, UWC, VWC, WWC, northdir);
		massFluxTC = WEqn->massFluxCorrection(PC, UWC, VWC, WWC, topdir);

		UC = UEqn->cellVelocityCorrection(PC, U_, V_, W_, eastdir);
		VC = VEqn->cellVelocityCorrection(PC, U_, V_, W_, northdir);
		WC = WEqn->cellVelocityCorrection(PC, U_, V_, W_, topdir);
		//********* end of correction of velocity and massflux******************//
//		cout<<" UWC "<<endl;
//		fieldOper.print3dmat(UWC);
//		cout<<" UWC end "<<endl;


//		cout<< " massFluxEC "<<endl;
//		fieldOper.print3dmat(massFluxEC);
//		cout<<" massFluxEC end " <<endl;	

		//*********** adding the correction********************//
		for( unsigned int i = 1; i<UW.size() - 1;i++)
		{
			for(unsigned int j = 1; j<UW[0].size() - 1; j++)
			{
				for(unsigned int k = 1 ; k< UW[0][0].size() - 1; k++)
				{
					UW[i][j][k].value += UWC[i][j][k].value;
					massFluxE[i][j][k].value += massFluxEC[i][j][k].value;
				}
			}
		}


		for( unsigned int i = 1; i<VW.size() -1; i++)
		{
			for(unsigned int j = 1; j<VW[0].size() -1; j++)
			{
				for(unsigned int k = 1; k< VW[0][0].size() -1; k++)
				{
					VW[i][j][k].value += VWC[i][j][k].value;
					massFluxN[i][j][k].value += massFluxNC[i][j][k].value;
				}
			}
		}


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


		for( unsigned int i = 1; i<P.size() - 1;i++)
		{
			for(unsigned int j = 1; j<P[0].size() - 1; j++)
			{
				for(unsigned int k = 1 ; k< P[0][0].size() - 1; k++)
				{
					P[i][j][k].value  = 0.8*P[i][j][k].value + 0.2*PC[i][j][k].value;
					U_[i][j][k].value = 0.2*U_[i][j][k].value + 0.8*UC[i][j][k].value;
					V_[i][j][k].value = 0.2*V_[i][j][k].value + 0.8*VC[i][j][k].value;
					W_[i][j][k].value = 0.2*W_[i][j][k].value + 0.8*WC[i][j][k].value;
				}
			}
		}

		forTopBoundary(W_)
		{
			W_[i][j][k].value = WW[i][j][k-1].value;

		}
		

		//****************extrapolate pressure**************************//
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, east);
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, west);
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, north);
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, south);
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, bottom);
		//***********end of pressure extrapolation***********************//


		U = U_;
		V = V_;
		W = W_;
		delete UEqn;
		delete VEqn;
		delete WEqn;
		delete PEqn;
		
//		cout<< "massFluxE at last"<<endl;
//		fieldOper.print3dmat(massFluxE);
//		cout<<" massFluxE end" <<endl;	
		cout<<"end of iteration "<<iter<<endl;

		fileWriter fileWriter_;
		string filename = "zrotatingmicrochannel";
		int time1 = iter;

		fileWriter_.writeUVWP(filename, time1, mymesh_,U,V,W,P);
	}
	

		//********* end of addition of correction***************//
		//*********************printing out the value *****************
	//	cout<<" U solved" <<endl;
	//	fieldOper.print3dmat(U);
		//cout<<" V partial solver" <<endl;
		//fieldOper.print3dmat(massFluxN);
		//cout<<" W partial sovled" <<endl;
		//fieldOper.print3dmat(massFluxT);

		//cout<<" size of UW "<<UW.size()<<" " <<UW[0].size()<<" "<<UW[0][0].size()<<endl;
		//cout<<NI<<" " <<NJ<<" "<<NK<<endl;

				
	
    return 0;

}
