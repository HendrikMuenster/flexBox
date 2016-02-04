#ifndef flexDualizedDataTerm_H
#define flexDualizedDataTerm_H

#include "flexTermDual.h"
#include "flexTildeMultiOperatorMultiDual.h"

template < typename T >
class flexDualizedDataTerm : public flexTildeMultiOperatorMultiDual<T>
{
private:


public: // basicDualizedDataterm(alpha,A,f,varargin)
	flexVector<T> f, fAlphaSigma;

	flexDualizedDataTerm(T _alpha, flexVector<flexMatrix<T>> _operatorList, flexVector<T> _f) : flexTildeMultiOperatorMultiDual(_alpha, 1)
	{
		operatorList.push_back(_operatorList[0]);
		operatorListT.push_back(_operatorList[0]);
		operatorListT[0].transpose();

		f = _f;
		fAlphaSigma = _f;

		mySigma.push_back(operatorList[0].getMaxRowSumAbs());
		myTau.push_back(operatorListT[0].getMaxRowSumAbs());
	};
};

#endif