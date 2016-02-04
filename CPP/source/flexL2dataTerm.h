#ifndef flexL2DataTerm_H
#define flexL2DataTerm_H

#include "flexTermPrimal.h"

template < typename T >
class flexL2DataTerm : public flexTermPrimal<T>
{
	private:
		flexVector<T> f, fAlphaTau;

	public:
		flexL2DataTerm(flexVector<int> _dims, T _alpha, flexVector<T> _f) : flexTermPrimal(1, _dims, _alpha)
		{
			f = _f;
			fAlphaTau = f;
		};

		void applyProx(flexBoxData<T> &data, flexVector<T> tau, flexVector<int> primalNumbers)
		{
			T factor = static_cast<T>(1) / (static_cast<T>(1) + alpha*tau[primalNumbers[0]]);

			fAlphaTau = f;
			fAlphaTau *= (alpha * tau[primalNumbers[0]]);

			data.x[primalNumbers[0]] = data.xTilde[primalNumbers[0]];
			data.x[primalNumbers[0]]+=(fAlphaTau);
			data.x[primalNumbers[0]] *= factor;

		};
};

#endif