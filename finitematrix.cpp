#include "finitematrix.hpp"
#include <iomanip> // std::setprecision

FiniteMatrix::FiniteMatrix():
value(0.0), apnot(0.0), ae(0.0), aw(0.0), an(0.0), as(0.0), at(0.0), ab(0.0), sp(0.0), sf(0.0) 
{

}




FiniteMatrix::~FiniteMatrix()
{

}



//operator overlaoding for +, - , *, &&

FiniteMatrix::finiteMat operator+(const FiniteMatrix::finiteMat& lhs, const FiniteMatrix::finiteMat& rhs)
{
	FiniteMatrix::finiteMat result(lhs);

	forAllInternal(lhs)
	{
		result[i][j][k].value += rhs[i][j][k].value;
		result[i][j][k].apnot += rhs[i][j][k].apnot;
		result[i][j][k].ae    += rhs[i][j][k].ae;
		result[i][j][k].aw    += rhs[i][j][k].aw;
		result[i][j][k].an    += rhs[i][j][k].an;
		result[i][j][k].as    += rhs[i][j][k].as;
		result[i][j][k].at    += rhs[i][j][k].at;
		result[i][j][k].ab    += rhs[i][j][k].ab;
		result[i][j][k].sp    += rhs[i][j][k].sp;
		result[i][j][k].sf    += rhs[i][j][k].sf;
	}
	return result;

}




// && element wise multiplication** *** NOT matrix multiplication *****
FiniteMatrix::finiteMat operator&&(const FiniteMatrix::finiteMat& lhs, const FiniteMatrix::finiteMat& rhs)
{
	FiniteMatrix::finiteMat result(lhs);

	forAllInternal(lhs)
	{
		result[i][j][k].value *= rhs[i][j][k].value;
		result[i][j][k].apnot *= rhs[i][j][k].apnot;
		result[i][j][k].ae    *= rhs[i][j][k].ae;
		result[i][j][k].aw    *= rhs[i][j][k].aw;
		result[i][j][k].an    *= rhs[i][j][k].an;
		result[i][j][k].as    *= rhs[i][j][k].as;
		result[i][j][k].at    *= rhs[i][j][k].at;
		result[i][j][k].ab    *= rhs[i][j][k].ab;
		result[i][j][k].sp    *= rhs[i][j][k].sp;
		result[i][j][k].sf    *= rhs[i][j][k].sf;

	}
	return result;

}





FiniteMatrix::finiteMat operator*(const double dblval, const FiniteMatrix::finiteMat& rhs)
{
	FiniteMatrix::finiteMat result(rhs);

	forAllInternal(result)
	{
		result[i][j][k].value *= dblval;
	}
	return result;

}




FiniteMatrix::finiteMat operator-(const FiniteMatrix::finiteMat& lhs, const FiniteMatrix::finiteMat& rhs)
{
	FiniteMatrix::finiteMat result(lhs);

	forAllInternal(lhs)
	{
		result[i][j][k].value -= rhs[i][j][k].value;
		result[i][j][k].apnot -= rhs[i][j][k].apnot;
		result[i][j][k].ae    -= rhs[i][j][k].ae;
		result[i][j][k].aw    -= rhs[i][j][k].aw;
		result[i][j][k].an    -= rhs[i][j][k].an;
		result[i][j][k].as    -= rhs[i][j][k].as;
		result[i][j][k].at    -= rhs[i][j][k].at;
		result[i][j][k].ab    -= rhs[i][j][k].ab;
		result[i][j][k].sf    -= rhs[i][j][k].sf;
		result[i][j][k].sp    -= rhs[i][j][k].sp;

	}
	return result;

}



//printing Ap =  value
void FiniteMatrix::print3dmat(finiteMat& vec)
{
	for(unsigned int i = 0; i<vec.size(); i++)
	{
		for(unsigned int j = 0 ; j<vec[0].size(); j++)
		{
			for(unsigned int k = 0; k< vec[0][0].size(); k++)
			{
				std::cout << vec[i][j][k].value<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
		cout<<endl;
	}
}




//printing source value
void FiniteMatrix::print3dsource(finiteMat& vec)
{
	for(unsigned int i = 0; i<vec.size(); i++)
	{
		for(unsigned int j = 0 ; j<vec[i].size(); j++)
		{
			for(unsigned int k = 0; k< vec[i][j].size(); k++)
			{
				std::cout <<std::setprecision(3)<< vec[i][j][k].sp<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
		cout<<endl;
	}
}
