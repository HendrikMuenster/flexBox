#ifndef flexL1dualizedOperatorAniso_H
#define flexL1dualizedOperatorAniso_H

#include "flexTermDual.h"
#include "flexTildeMultiOperatorMultiDual.h"

template < typename T >
class flexL1dualizedOperatorAniso : public flexTildeMultiOperatorMultiDual<T>
{
private:

public:
	flexL1dualizedOperatorAniso(T _alpha, flexVector<flexMatrix<T>> _operatorList) : flexTildeMultiOperatorMultiDual(_alpha, _operatorList.size())
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
		for (int j = 0; j < operatorList.size(); ++j)
		{
			data.yTilde[dualNumbers[j]].project_minMax(alpha);
			data.y[dualNumbers[j]] = data.yTilde[dualNumbers[j]];
		}
	};
};

#endif