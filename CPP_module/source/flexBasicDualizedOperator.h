#ifndef flexBasicdualizedOperator_H
#define flexBasicdualizedOperator_H

#include "flexTermDual.h"
#include "flexTildeMultiOperatorMultiDual.h"

template < typename T >
class flexBasicDualizedOperator : public flexTildeMultiOperatorMultiDual<T>
{
private:

public:
	flexBasicDualizedOperator(T _alpha, int numberPrimals, flexVector<flexMatrix<T> > _operatorList) : flexTildeMultiOperatorMultiDual<T>(_alpha, numberPrimals, _operatorList.size() / numberPrimals)
	{
			
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