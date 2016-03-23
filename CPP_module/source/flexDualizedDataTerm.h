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

	flexDualizedDataTerm(T _alpha, flexVector<flexMatrix<T> > _operatorList, flexVector<T> _f) : flexTildeMultiOperatorMultiDual<T>(_alpha, 1, 1)
	{
		f = _f;
		fAlphaSigma = _f;

		int numberPrimals = _operatorList.size();

		this->operatorList = _operatorList;

		//create sigma and tau
		this->myTau.resize(numberPrimals, 0.0);
		this->mySigma.resize(_operatorList.size() / numberPrimals, 0.0);

		for (int i = 0; i < _operatorList.size() / numberPrimals; ++i)
		{
			for (int j = 0; j < numberPrimals; ++j)
			{
				int opNum = i*numberPrimals + j;

				this->operatorListT.push_back(this->operatorList[opNum]);
				this->operatorListT[opNum].transpose();

				this->mySigma[i] += this->operatorList[opNum].getMaxRowSumAbs();
				this->myTau[j] += this->operatorListT[opNum].getMaxRowSumAbs();
			}
		}




	};
};

#endif