#ifndef flexL1dualizedOperatorAniso_H
#define flexL1dualizedOperatorAniso_H

#include "flexTermDual.h"
#include "flexBasicDualizedOperator.h"

template < typename T >
class flexL1dualizedOperatorAniso : public flexBasicDualizedOperator<T>
{
	public:
		flexL1dualizedOperatorAniso(T _alpha, int numberPrimals, flexVector<flexLinearOperator<T>* > _operatorList) : flexBasicDualizedOperator<T>(_alpha, numberPrimals, _operatorList){};

		void applyProx(flexBoxData<T> &data, const flexVector<T> &sigma, const flexVector<int> &dualNumbers, const flexVector<int> &primalNumbers)
		{
			for (int j = 0; j < dualNumbers.size(); ++j)
			{
				data.yTilde[dualNumbers[j]].project_minMax(this->alpha);
				data.y[dualNumbers[j]] = data.yTilde[dualNumbers[j]];
			}
		};
};

#endif