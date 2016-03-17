#ifndef flexFrobeniusDualizedOperator_H
#define flexFrobeniusDualizedOperator_H

#include "flexTermDual.h"
#include "flexBasicDualizedOperator.h"

template < typename T >
class flexFrobeniusDualizedOperator : public flexBasicDualizedOperator<T>
{
public:
	flexFrobeniusDualizedOperator(T _alpha, int numberPrimals, flexVector<flexMatrix<T>> _operatorList) : flexBasicDualizedOperator(_alpha, numberPrimals, _operatorList) {};

	void applyProx(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		int numElements = data.y[dualNumbers[0]].size();

		flexVector<T> normTmp(numElements, static_cast<T>(0));
		T normTmp2 = static_cast<T>(0);

		//for all dual variables
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			normTmp = data.yTilde[dualNumbers[i]];
			normTmp.pow2();

			normTmp2 += normTmp.sum();
		}

		normTmp2 = max(static_cast<T>(1), std::sqrt(normTmp2) / alpha);

		for (int i = 0; i < numElements; ++i)
		{
			for (int j = 0; j < dualNumbers.size(); ++j)
			{
				data.y[dualNumbers[j]][i] = data.yTilde[dualNumbers[j]][i] / normTmp2;
			}
		}
	};
};

#endif