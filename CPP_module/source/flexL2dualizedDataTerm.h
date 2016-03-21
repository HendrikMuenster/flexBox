#ifndef flexL2DualizedDataTerm_H
#define flexL2DualizedDataTerm_H

#include "flexTermDual.h"
#include "flexDualizedDataTerm.h"

template < typename T >
class flexL2dualizedDataTerm : public flexDualizedDataTerm<T>
{
private:
	bool initiated;
	T factor;

public:
	flexL2dualizedDataTerm(T _alpha, flexVector<flexMatrix<T>> _operatorList, flexVector<T> _f) : flexDualizedDataTerm(_alpha, _operatorList, _f) 
	{
		initiated = false;
	};

	void initiate(T sigma)
	{
		factor = alpha / (alpha + sigma); 
		
		fAlphaSigma = f;
		fAlphaSigma *= sigma;

		initiated = true;
	}

	void applyProx(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		//main.y{dualNumber} = (obj.factor/(main.params.sigma{dualNumber}+obj.factor)) * (main.yTilde{dualNumber} - main.params.sigma{dualNumber}*obj.f(:));

		if (!initiated)
		{
			initiate(sigma[dualNumbers[0]]);
		}
		
		data.y[dualNumbers[0]] = data.yTilde[dualNumbers[0]];
		data.y[dualNumbers[0]] -= fAlphaSigma;
		data.y[dualNumbers[0]] *= factor;
	};
};

#endif