#ifndef flexL1DualizedOperator_H
#define flexL1DualizedOperator_H

#include "flexTermDual.h"
#include "flexDualizedOperator.h"

template < typename T >
class flexL1DualizedOperator : public flexDualizedOperator<T>
{
private:

public:
	flexL1DualizedDataterm(T _alpha, flexMatrix<T> _operator, T _myTau, T _mySigma) : flexDualizedOperator(_alpha, _operator, _myTau, _mySigma)
	{

	};

	void applyProx(flexBoxData<T> &data, T sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		// main.y{dualNumber} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber}));

		data.y[dualNumbers[0]] = data.yTilde[dualNumbers[0]];
		data.y[dualNumbers[0]].project_minMax();
	};
};

#endif