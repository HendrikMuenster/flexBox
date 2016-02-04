#ifndef flexL2dualizedOperator_H
#define flexL2dualizedOperator_H

#include "flexTermDual.h"
#include "flexTildeMultiOperatorMultiDual.h"

template < typename T >
class flexL2dualizedOperator : public flexTildeMultiOperatorMultiDual<T>
{
	private:

	public:

	flexL2dualizedOperator(T _alpha, flexVector<flexMatrix<T>> _operatorList) : flexTildeMultiOperatorMultiDual(_alpha, _operatorList.size())
	{
		operatorList = _operatorList;

		for (int i = 0; i < operatorList.size(); ++i)
		{
			operatorListT.push_back(operatorList[i]);
			operatorListT[i].transpose();

			mySigma.push_back(operatorList[i].getMaxRowSumAbs());
			myTau.push_back(operatorListT[i].getMaxRowSumAbs());
		}
	};

	void applyProx(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{


		for (int i = 0; i < operatorList.size(); ++i)
		{
			T factor = alpha / (alpha + sigma[dualNumbers[i]]);

			// should be data.y[dualNumbers[i]] = data.yTilde[dualNumbers[i]] * factor;
			data.y[dualNumbers[i]] = data.yTilde[dualNumbers[i]];
			data.y[dualNumbers[i]] *= factor;
		}
		
	};
};

#endif