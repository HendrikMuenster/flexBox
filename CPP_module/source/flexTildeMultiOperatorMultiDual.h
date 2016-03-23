#ifndef flexTildeMultiOperatorMultiDual_H
#define flexTildeMultiOperatorMultiDual_H

#include "flexVector.h"
#include "flexBox.h"

template < typename T >
class flexTildeMultiOperatorMultiDual : public flexTermDual<T>
{
private:


public:
	flexTildeMultiOperatorMultiDual(T _alpha, int numberPrimals, int _numberVars) : flexTermDual<T>(_alpha, numberPrimals, _numberVars){};

	void yTilde(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		//set yTilde = y;
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			data.yTilde[dualNumbers[i]] = data.y[dualNumbers[i]];
		}

		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				//Check here
				int operatorNumber = primalNumbers.size() * i + j;

				flexVector<T> xBarTmp(data.xBar[primalNumbers[j]]);
				xBarTmp *= sigma[dualNumbers[i]];

				this->operatorList[operatorNumber].timesPlus(xBarTmp,data.yTilde[dualNumbers[i]]);

				//data.xBar[primalNumbers[j]].print();
				//xBarTmp.print();
				//data.yTilde[dualNumbers[i]].print();
			}
		}
	}

	void xTilde(flexBoxData<T> &data, flexVector<T> tau, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				//Check here
				int operatorNumber = primalNumbers.size() * i + j;

				flexVector<T> yTmp(data.y[dualNumbers[i]]);
				yTmp *= tau[primalNumbers[j]];

				this->operatorListT[operatorNumber].timesMinus(yTmp,data.xTilde[primalNumbers[j]]);
			}
		}
	}

	void yError(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				//Check here
				int operatorNumber = primalNumbers.size() * i + j;

				flexVector<T> xTmp(data.x[primalNumbers[j]]);
				xTmp -= data.xOld[primalNumbers[j]];

				this->operatorList[operatorNumber].timesPlus(xTmp, data.yError[dualNumbers[i]]);
			}
		}
	}

	void xError(flexBoxData<T> &data, flexVector<T> tau, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				//Check here
				int operatorNumber = primalNumbers.size() * i + j;

				flexVector<T> yTmp(data.y[dualNumbers[i]]);
				yTmp -= data.yOld[dualNumbers[i]];

				this->operatorListT[operatorNumber].timesMinus(yTmp, data.xError[primalNumbers[j]]);
			}
		}
	}
};

#endif