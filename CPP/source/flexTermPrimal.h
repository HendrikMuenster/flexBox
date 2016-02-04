#ifndef flexTermPrimal_H
#define flexTermPrimal_H

#include "flexBoxData.h"
#include "flexVector.h"

template < typename T > class flexBoxData;

template < typename T >
class flexTermPrimal
{
	private:
		int numberVars;
		flexVector<int> dims;
		
	public:
		T alpha;

		flexTermPrimal(int _numberVars, flexVector<int> _dims, T _alpha)
		{
			numberVars = _numberVars;
			dims = _dims;
			alpha = _alpha;
		};

		int getNumberVars()
		{
			return numberVars;
		}

		flexVector<int> getDims()
		{
			return dims;
		}

		virtual void applyProx(flexBoxData<T> &data, flexVector<T> tau, flexVector<int> primalNumbers) = 0;
};

#endif