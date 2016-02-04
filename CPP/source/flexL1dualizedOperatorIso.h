#ifndef flexL1dualizedOperatorIso_H
#define flexL1dualizedOperatorIso_H

#include "flexTermDual.h"
#include "flexTildeMultiOperatorMultiDual.h"

template < typename T >
class flexL1dualizedOperatorIso : public flexTildeMultiOperatorMultiDual<T>
{
private:

public:
	flexL1dualizedOperatorIso(T _alpha, flexVector<flexMatrix<T>> _operatorList) : flexTildeMultiOperatorMultiDual(_alpha, _operatorList.size())
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
		flexVector<T> normTmp(data.y[dualNumbers[0]].size(), static_cast<T>(0));
		flexVector<T> normTmp2(data.y[dualNumbers[0]].size(), static_cast<T>(0));

		for (int i = 0; i < operatorList.size(); ++i)
		{
			normTmp = data.yTilde[dualNumbers[i]];
			normTmp *= normTmp;

			normTmp2 += normTmp;
		}

		for (int i = 0; i < data.y[dualNumbers[0]].size(); ++i)
		{
			normTmp2[i] = max(static_cast<T>(1), sqrt(normTmp2[i]) / alpha);

			for (int j = 0; j < operatorList.size(); ++j)
			{
				data.y[dualNumbers[j]][i] = data.yTilde[dualNumbers[j]][i] / normTmp2[i];
			}
		}

	};
};

#endif