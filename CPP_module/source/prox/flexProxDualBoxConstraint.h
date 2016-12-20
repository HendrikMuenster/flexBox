#ifndef flexProxDualBoxConstraint_H
#define flexProxDualBoxConstraint_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualBoxConstraint : public flexProx<T, Tvector>
{
private:
	T minVal;
	T maxVal;
public:

	flexProxDualBoxConstraint(T _minVal, T _maxVal) : flexProx<T, Tvector>(dualBoxConstraintProx)
	{
		minVal = _minVal;
		maxVal = _maxVal;
	}

	~flexProxDualBoxConstraint()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		#if __CUDACC__
		#else
			for (int k = 0; k < dualNumbers.size(); k++)
			{
				T* ptrY = data->y[dualNumbers[k]].data();
				T* ptrYtilde = data->yTilde[dualNumbers[k]].data();
				T* ptrSigma = data->sigmaElt[dualNumbers[k]].data();

				int numElements = (int)data->yTilde[dualNumbers[k]].size();

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					ptrY[i] = myMax<T>((T)0, ptrYtilde[i] - ptrSigma[i] * this->maxVal) + myMin<T>((T)0, ptrYtilde[i] - ptrSigma[i] * this->minVal);
				}
			}
		#endif
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<std::vector<T>> &fList)
	{

	}
};

#endif