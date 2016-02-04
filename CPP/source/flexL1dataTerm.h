#ifndef flexL1DataTerm_H
#define flexL1DataTerm_H

#include "flexTermPrimal.h"

template < typename T >
class flexL1DataTerm : public flexTermPrimal<T>
{
private:
	flexVector<T> f, fAlphaTau;

public:
	flexL1DataTerm(flexVector<int> _dims, T _alpha, flexVector<T> _f) : flexTermPrimal(1, _dims, _alpha)
	{
		f = _f;
	};

	void applyProx(flexBoxData<T> &data, flexVector<T> tau, flexVector<int> primalNumbers)
	{
		T factor = alpha*tau[primalNumbers[0]];

		//T *varPointer = &data.xTilde[primalNumbers[0]];

		for (int i = 0; i < data.xTilde[primalNumbers[0]].size(); ++i)
		{
			if ((data.xTilde[primalNumbers[0]][i] - f[i]) < -factor)
			{
				data.x[primalNumbers[0]][i] = data.xTilde[primalNumbers[0]][i] + factor;
			}
			else if ((data.xTilde[primalNumbers[0]][i] - f[i]) > factor)
			{
				data.x[primalNumbers[0]][i] = data.xTilde[primalNumbers[0]][i] - factor;
			}
			else
			{
				data.x[primalNumbers[0]][i] = f[i];
			}
		}
	};
};

#endif