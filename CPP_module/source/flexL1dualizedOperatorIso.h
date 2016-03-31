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
		flexL1dualizedOperatorIso(T _alpha, int numberPrimals, flexVector<flexLinearOperator<T>* > _operatorList) : flexBasicDualizedOperator<T>(_alpha, numberPrimals, _operatorList)
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

			int elements = normTmp.size();

			//for all dual variables
			for (int i = 0; i < dualNumbers.size(); ++i)
			{
				//printf("Dual %d norm\n", dualNumbers[i]);

				normTmp = data.yTilde[dualNumbers[i]];
				normTmp.pow2();
				normTmp2 += normTmp;
			}

			T numberOne = static_cast<T>(1);

			#pragma omp parallel for
			for (int j = 0; j < elements; ++j)
			{
				normTmp2[j] = myMax(numberOne, sqrt(normTmp2[j]) / this->alpha);
			}

			for (int i = 0; i < dualNumbers.size(); ++i)
			{
				//printf("Dual %d proj\n", dualNumbers[i]);

				data.y[dualNumbers[i]] = data.yTilde[dualNumbers[i]];

				data.y[dualNumbers[i]].vectorDivide(normTmp2);
			}

		};
};

#endif