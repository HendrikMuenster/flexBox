#ifndef flexL1dualizedOperatorIso_H
#define flexL1dualizedOperatorIso_H

#include "flexTermDual.h"
#include "flexBasicDualizedOperator.h"

template < typename T >
class flexL1dualizedOperatorIso : public flexBasicDualizedOperator<T>
{
	private:
		flexVector<T> normTmp, normTmp2;
		int numElements;
		bool initiated;

	public:
		flexL1dualizedOperatorIso(T _alpha, int numberPrimals, flexVector<flexMatrix<T>> _operatorList) : flexBasicDualizedOperator(_alpha, numberPrimals, _operatorList)
		{
			initiated = false;
		};

		void initiate(int elements)
		{
			int numElements = elements;

			normTmp.resize(numElements, static_cast<T>(0));
			normTmp2.resize(numElements, static_cast<T>(0));

			initiated = true;
		}

		void applyProx(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
		{
			if (!initiated)
			{
				initiate(data.y[dualNumbers[0]].size());
			}
			
			normTmp2.scalarEquals(static_cast<T>(0));

			//for all dual variables
			for (int i = 0; i < dualNumbers.size(); ++i)
			{
				normTmp = data.yTilde[dualNumbers[i]];
				normTmp.pow2();

				normTmp2 += normTmp;
			}

			normTmp2.sqrt();
			normTmp2.scalarDivide(alpha);
			normTmp2.scalarMax(static_cast<T>(1));

			for (int j = 0; j < dualNumbers.size(); ++j)
			{
				data.y[dualNumbers[j]] = data.yTilde[dualNumbers[j]];

				data.y[dualNumbers[j]].vectorDivide(normTmp2);
			}

		};
};

#endif