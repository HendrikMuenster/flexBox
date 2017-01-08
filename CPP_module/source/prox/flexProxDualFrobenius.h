#ifndef flexProxDualFrobenius_H
#define flexProxDualFrobenius_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualFrobenius : public flexProx<T, Tvector>
{
private:

public:

	flexProxDualFrobenius() : flexProx<T, Tvector>(dualFrobeniusProx)
	{
	}

	~flexProxDualFrobenius()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		#ifdef __CUDACC__
            printf("flexProxDualFrobenius Prox not implemented for CUDA\n");
		#else
			T norm = (T)0;
			for (int k = 0; k < dualNumbers.size(); k++)
			{
				T* ptrYTilde = data->yTilde[dualNumbers[k]].data();

				int numElements = (int)data->yTilde[dualNumbers[k]].size();

				#pragma omp parallel for reduction(+: norm)
				for (int i = 0; i < numElements; i++)
				{
					norm += ptrYTilde[i] * ptrYTilde[i];
				}
			}

			norm = (T)1 / std::max((T)1, std::sqrt(norm) / alpha);

			for (int k = 0; k < dualNumbers.size(); k++)
			{
				T* ptrY = data->y[dualNumbers[k]].data();
				T* ptrYTilde = data->yTilde[dualNumbers[k]].data();

				int numElements = (int)data->yTilde[dualNumbers[k]].size();

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					ptrY[i] = ptrYTilde[i] * norm;
				}
			}
		#endif
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{

	}
};

#endif