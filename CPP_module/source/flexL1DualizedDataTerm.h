#ifndef flexL1DualizedDataTerm_H
#define flexL1DualizedDataTerm_H

#include "flexTermDual.h"
#include "flexDualizedDataTerm.h"

template < typename T >
class flexL1dualizedDataTerm : public flexDualizedDataTerm<T>
{
private:
	bool initiated;
	
public:
	flexL1dualizedDataTerm(T _alpha, flexVector<flexLinearOperator<T>* > _operatorList, flexVector<T> _f) : flexDualizedDataTerm<T>(_alpha, _operatorList, _f)
	{
		initiated = false;
	};

	void initiate(T sigma)
	{
		this->fAlphaSigma = this->f;
		this->fAlphaSigma.scalarMult(sigma*this->alpha);

		initiated = true;
	}

	void applyProx(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		// main.y{dualNumber} = max(-obj.factor,min(obj.factor,main.yTilde{dualNumber} - main.params.sigma{dualNumber}*obj.f(:)));

		if (!initiated)
		{
			initiate(sigma[dualNumbers[0]]);
		}

		//printf("Dual %d direct\n", dualNumbers[0]);

		data.y[dualNumbers[0]] = data.yTilde[dualNumbers[0]];

		data.y[dualNumbers[0]] -= this->fAlphaSigma;

		data.y[dualNumbers[0]].project_minMax(this->alpha);
	};
};

#endif