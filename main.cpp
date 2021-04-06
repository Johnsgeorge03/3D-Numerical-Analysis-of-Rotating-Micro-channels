#include<iostream>
#include "mesh.hpp"
#include "solution.hpp"
#include "foralloperations.hpp"
#include "fields.hpp"
#include "finitematrix.hpp"
#include "equation.hpp"
#include "finitevolumeoperations.hpp"
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
	#include "initializevariables.hpp"
    	#include "initializefinitematrixvar.hpp" 

	int iter = 0;
	int maxit = 20;
	for(iter = 0; iter< maxit; iter++)
	{
		//update pressure boundary condition
		//*******left out for now**************
		

		int eastdir = 1;
		int northdir = 2; 
		int topdir = 3;
		unsigned int iterations = 50;


		
		//************Coefficient calculation**********************//
		Equation *UEqn = new Equation
		(
			fvm::convDiffusiveTerm(U, massFluxE, massFluxN, massFluxT) 
			+ fvm::pressureGrad(P, eastdir)
			+ fvm::forceSource(W, sol, eastdir)
			+ fvm::timeCoeff(U, sol)

		);
		
		UEqn->assembleEquation(massFluxE, massFluxN, massFluxT);
		
		Equation *VEqn = new Equation
		(
			fvm::convDiffusiveTerm(V, massFluxE, massFluxN, massFluxT)
			+ fvm::pressureGrad(P, northdir)
			+ fvm::forceSource(U, sol, northdir) // U is not used in calc
			+ fvm::timeCoeff(V, sol)
		);

		VEqn->assembleEquation(massFluxE, massFluxN, massFluxT);

		Equation *WEqn = new Equation
		(
			fvm::convDiffusiveTerm(W, massFluxE, massFluxN, massFluxT)
			+ fvm::pressureGrad(P, topdir)
			+ fvm::forceSource(U, sol, topdir)
			+ fvm::timeCoeff(W, sol)
		);

		WEqn->assembleEquation(massFluxE, massFluxN, massFluxT);
		//**************end of coefficient calculation********************//



		//**************Solving the equation that satisfy only momemtum equation**************//
		Fields::vec3dField U_ = UEqn->solveVelocity(U, iterations);
		Fields::vec3dField V_ = VEqn->solveVelocity(V, iterations);
		Fields::vec3dField W_ = WEqn->solveVelocity(W, iterations);
		fieldOper.copyOutletVelocity(U_);
		fieldOper.copyOutletVelocity(V_);
		fieldOper.copyOutletVelocity(W_);
		//***********************end of solving partial equation******************************//
		


		//****************momentum interpolation****************************//
		UW = UEqn->momentumInterpolation(U_, P, eastdir);
		VW = VEqn->momentumInterpolation(V_, P, northdir);
		WW = WEqn->momentumInterpolation(W_, P, topdir);
		//***************end of momentum interpolation**********************//	
		



		//****************mass flux calculation*****************************//
		// additional vel field is sent to extract the density and face area//
		massFluxE = UEqn->massFluxCalculation(UW, U, eastdir);
		massFluxN = VEqn->massFluxCalculation(VW, V, northdir);
		massFluxT = WEqn->massFluxCalculation(WW, W, topdir);
		//**************end of mass flux calculation***********************//





		//***************pressure coefficient assembly*********************//
		Equation *PEqn =  new Equation
		(
			fvm::pressureCorrectionCoeff(P, massFluxE, massFluxN, massFluxT)
			+ fvm::pressureCorrectionSource(P, massFluxE, massFluxN, massFluxT)
		);

		PEqn->assemblePressureEquation();
		//**************end of pressure coefficient assembly*******************//




		//*****************pressure correction solve************************//
		PC =  PEqn->solvePressure(P, iterations);
		//*************end of pressure correction solve*********************//




		//****************extrapolate pressure**************************//
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, east);
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, west);
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, north);
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, south);
		fieldOper.linearextrapolateCondition(P, mymesh_.FX, mymesh_.FY, mymesh_.FZ, bottom);
		//***********end of pressure extrapolation***********************//





		//*******************correction of velocity and massflux************//
		UWC = UEqn->velocityCorrection(PC, UW, VW, WW, eastdir);
		VWC = VEqn->velocityCorrection(PC, UW, VW, WW, northdir);
		WWC = WEqn->velocityCorrection(PC, UW, VW, WW, topdir);

		massFluxEC = UEqn->massFluxCorrection(PC, UWC, VWC, WWC, eastdir);
		massFluxNC = VEqn->massFluxCorrection(PC, UWC, VWC, WWC, northdir);
		massFluxTC = WEqn->massFluxCorrection(PC, UWC, VWC, WWC, topdir);
		//********* end of correction of velocity and massflux******************//



		//*********** adding the correction********************//
		for( unsigned int i = 0; i<UW.size();i++)
		{
			for(unsigned int j = 0; j<UW[0].size(); j++)
			{
				for(unsigned int k = 0 ; k< UW[0][0].size(); k++)
				{
					UW[i][j][k].value += UWC[i][j][k].value;
					massFluxE[i][j][k].value += massFluxEC[i][j][k].value;
				}
			}
		}

		for( unsigned int i = 0; i<VW.size();i++)
		{
			for(unsigned int j = 0; j<VW[0].size(); j++)
			{
				for(unsigned int k = 0 ; k< VW[0][0].size(); k++)
				{
					VW[i][j][k].value += VWC[i][j][k].value;
					massFluxN[i][j][k].value += massFluxNC[i][j][k].value;
				}
			}
		}

		for( unsigned int i = 0; i<WW.size();i++)
		{
			for(unsigned int j = 0; j<WW[0].size(); j++)
			{
				for(unsigned int k = 0 ; k< WW[0][0].size(); k++)
				{
					WW[i][j][k].value += WWC[i][j][k].value;
					massFluxT[i][j][k].value += massFluxTC[i][j][k].value;
				}
			}
		}

		for( unsigned int i = 0; i<P.size();i++)
		{
			for(unsigned int j = 0; j<P[0].size(); j++)
			{
				for(unsigned int k = 0 ; k< P[0][0].size(); k++)
				{
					P[i][j][k].value += 0.8*PC[i][j][k].value;
				}
			}
		}

		cout<<"end of iteration "<<iter<<endl;
	}
	

		//********* end of addition of correction***************//
		//*********************printing out the value *****************
		//cout<<" U solved" <<endl;
		//fieldOper.print3dmat(massFluxE);
		//cout<<" V partial solver" <<endl;
		//fieldOper.print3dmat(massFluxN);
		//cout<<" W partial sovled" <<endl;
		//fieldOper.print3dmat(massFluxT);

		//cout<<" size of UW "<<UW.size()<<" " <<UW[0].size()<<" "<<UW[0][0].size()<<endl;
		//cout<<NI<<" " <<NJ<<" "<<NK<<endl;

		cout<< " V "<< endl;
		for(unsigned int i = 0; i<U.size(); i++)
		{
			for(unsigned int j = 0 ; j<U[i].size(); j++)
			{
				for(unsigned int k = 0; k< U[i][j].size(); k++)
				{
					std::cout << V[i][j][k].value<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
			cout<<endl;
		}
		
	
    return 0;

}
