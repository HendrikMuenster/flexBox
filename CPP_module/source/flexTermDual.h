#ifndef flexTermDual_H
#define flexTermDual_H

#include "flexBoxData.h"
#include "flexVector.h"
#include "flexLinearOperator.h"

//template < typename T > class flexBoxData;

template < typename T >
class flexTermDual
{
	private:
		int numberVars;
		int numberPrimals;

	public:
		T alpha;
		flexVector<T> myTau;
		flexVector<T> mySigma;
		flexVector<flexLinearOperator<T>* > operatorList;
		flexVector<flexLinearOperator<T>* > operatorListT;

		flexTermDual(T _alpha, int _numberPrimals, int _numberVars)
		{
			numberVars = _numberVars;
			numberPrimals = _numberPrimals;
			alpha = _alpha;
		}

		int getNumberVars()
		{
			return numberVars;
		}

		int dualVarLength(int num)
		{
			return operatorList[num]->getNumRows();
		}
		
		virtual void applyProx(flexBoxData<T> &data, const flexVector<T> &sigma, const flexVector<int> &dualNumbers, const flexVector<int> &primalNumbers) = 0;

		virtual void yTilde(flexBoxData<T> &data, const flexVector<T> &sigma, const flexVector<int> &dualNumbers, const flexVector<int> &primalNumbers) = 0;

		virtual void xTilde(flexBoxData<T> &data, const flexVector<T> &tau, const flexVector<int> &dualNumbers, const flexVector<int> &primalNumbers) = 0;

		virtual void yError(flexBoxData<T> &data, const flexVector<T> &sigma, const flexVector<int> &dualNumbers, const flexVector<int> &primalNumbers) = 0;

		virtual void xError(flexBoxData<T> &data, const flexVector<T> &tau, const flexVector<int> &dualNumbers, const flexVector<int> &primalNumbers) = 0;
};

#endif