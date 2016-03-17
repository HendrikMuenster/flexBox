#ifndef flexKLDualizedDataTerm_H
#define flexKLDualizedDataTerm_H

#include "flexTermDual.h"
#include "flexDualizedDataTerm.h"

template < typename T >
class flexKLdualizedDataTerm : public flexDualizedDataTerm<T>
{
private:

public:
	flexKLdualizedDataTerm(T _alpha, flexVector<flexMatrix<T>> _operatorList, flexVector<T> _f) : flexDualizedDataTerm(_alpha, _operatorList, _f) {};

	void applyProx(flexBoxData<T> &data, flexVector<T> sigma, flexVector<int> dualNumbers, flexVector<int> primalNumbers)
	{
		//implements MATLAB code:
		//main.y{dualNumber} = 0.5*(1 + main.yTilde{dualNumber} - sqrt( (main.yTilde{dualNumber}-1).^2 + 4*main.params.sigma{dualNumber}*obj.f(:) ));


		flexVector<T> tmpTilde = data.yTilde[dualNumbers[0]];
		tmpTilde -= static_cast<T>(1);
		tmpTilde.pow2();

		fAlphaSigma = f;
		fAlphaSigma *= sigma[dualNumbers[0]];
		fAlphaSigma *= static_cast<T>(4);

		tmpTilde += fAlphaSigma;
		tmpTilde.sqrt();

		data.y[dualNumbers[0]] = data.yTilde[dualNumbers[0]];
		data.y[dualNumbers[0]] += static_cast<T>(1);
		data.y[dualNumbers[0]] -= tmpTilde;
		data.y[dualNumbers[0]] *= static_cast<T>(0.5);




	};
};

#endif