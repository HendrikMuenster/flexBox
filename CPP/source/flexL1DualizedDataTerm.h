#ifndef flexL1DualizedDataTerm_H
#define flexL1DualizedDataTerm_H

#include "flexTermDual.h"
#include "flexDualizedDataTerm.h"

template < typename T >
class flexL1dualizedDataTerm : public flexDualizedDataTerm<T>
{
private:
	
public:
	flexL1dualizedDataTerm(T _alpha, flexVector<flexMatrix<T>> _operatorList, flexVector<T> _f) : flexDualizedDataTerm(_alpha, _operatorList, _f)
	{

	};

	void applyProx(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		// main.y{dualNumber} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber} - main.params.sigma{dualNumber}*obj.f(:)));

		data.y[dualNumbers[0]] = data.yTilde[dualNumbers[0]];

		fAlphaSigma = f;
		fAlphaSigma *= alpha*sigma[dualNumbers[0]];
		data.y[dualNumbers[0]] -= fAlphaSigma;

		data.y[dualNumbers[0]].project_minMax(alpha);
	};
};

#endif