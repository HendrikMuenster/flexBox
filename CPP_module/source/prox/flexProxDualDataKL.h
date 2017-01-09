#ifndef flexProxDualKL_H
#define flexProxDualKL_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualDataKL : public flexProx<T, Tvector>
{
private:

public:

	flexProxDualDataKL() : flexProx<T, Tvector>(dualKLDataProx)
	{
	}

	~flexProxDualDataKL()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{
		#ifdef __CUDACC__
            printf("flexProxDualDataKL Prox not implemented for CUDA\n");
		#else
			for (int i = 0; i < dualNumbers.size(); i++)
			{
				T* ptrY = data->y[dualNumbers[i]].data();
				T* ptrYtilde = data->yTilde[dualNumbers[i]].data();
				T* ptrSigma = data->sigmaElt[dualNumbers[i]].data();

				T* ptrF = fList[i].data();

				int numElements = (int)data->yTilde[dualNumbers[i]].size();

				#pragma omp parallel for
				for (int j = 0; j < numElements; j++)
				{
					ptrY[j] = (T)0.5 * ((T)1 + ptrYtilde[j] - std::sqrt(myPow2<T>(ptrYtilde[j] - (T)1) + (T)4 * ptrSigma[j] * ptrF[j]));
				}
			}
		#endif
	}
};

#endif