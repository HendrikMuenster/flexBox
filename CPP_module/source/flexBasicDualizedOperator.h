#ifndef flexBasicdualizedOperator_H
#define flexBasicdualizedOperator_H

#include "flexTermDual.h"
#include "flexTildeMultiOperatorMultiDual.h"

template < typename T >
class flexBasicDualizedOperator : public flexTildeMultiOperatorMultiDual<T>
{
private:

public:
	flexBasicDualizedOperator(T _alpha, int numberPrimals, flexVector<flexMatrix<T>> _operatorList) : flexTildeMultiOperatorMultiDual(_alpha, numberPrimals, _operatorList.size() / numberPrimals)
	{
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