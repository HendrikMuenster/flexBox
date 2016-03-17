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

	flexDualizedDataTerm(T _alpha, flexVector<flexMatrix<T>> _operatorList, flexVector<T> _f) : flexTildeMultiOperatorMultiDual(_alpha, 1, 1)
	{
		f = _f;
		fAlphaSigma = _f;

		int numberPrimals = _operatorList.size();

		operatorList = _operatorList;

		//create sigma and tau
		myTau.resize(numberPrimals, 0.0);
		mySigma.resize(_operatorList.size() / numberPrimals, 0.0);

		for (int i = 0; i < _operatorList.size() / numberPrimals; ++i)
		{
			for (int j = 0; j < numberPrimals; ++j)
			{
				int opNum = i*numberPrimals + j;

				operatorListT.push_back(operatorList[opNum]);
				operatorListT[opNum].transpose();

				mySigma[i] += operatorList[opNum].getMaxRowSumAbs();
				myTau[j] += operatorListT[opNum].getMaxRowSumAbs();
			}
		}




	};
};

#endif