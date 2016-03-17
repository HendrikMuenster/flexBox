#ifndef flexDualizedOperator_H
#define flexDualizedOperator_H

#include "flexTermDual.h"
#include "flexTildeMultiOperatorMultiDual.h"

template < typename T >
class flexDualizedOperator : public flexTildeMultiOperatorMultiDual<T>
{
private:


public:

	flexDualizedOperator(T _alpha, flexMatrix<T> _operator, T _myTau, T _mySigma) : flexTildeMultiOperatorMultiDual(_alpha, 1)
	{
		operatorList.pushBack(_operator);
		operatorListT.pushBack(_operator);
		operatorListT[0].transpose();

		myTau = _myTau;
		mySigma = _mySigma;
	};
};

#endif     