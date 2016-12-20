#ifndef flexProxL2_H
#define flexProxL2_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualL2 : public flexProx<T, Tvector>
{
private:

public:

	flexProxDualL2() : flexProx<T, Tvector>(dualL2Prox)
	{
	}

	~flexProxDualL2()
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

				T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

				int numElements = (int)data->yTilde[dualNumbers[k]].size();

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					ptrY[i] = alpha / (ptrSigma[i] + alpha) * ptrYtilde[i];
				}
			}
		#endif
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<std::vector<T>> &fList)
	{

	}
};

#endif