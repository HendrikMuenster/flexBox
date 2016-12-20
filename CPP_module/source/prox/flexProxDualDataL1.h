#ifndef flexProxDualL1_H
#define flexProxDualL1_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualDataL1 : public flexProx<T, Tvector>
{
private:

public:

	flexProxDualDataL1() : flexProx<T, Tvector>(dualL1DataProx)
	{
	}

	~flexProxDualDataL1()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<std::vector<T>> &fList)
	{
		#if __CUDACC__
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
					ptrY[j] = myMin<T>(alpha, myMax<T>(-alpha, ptrYtilde[j] - alpha * ptrSigma[j] * ptrF[j]));
				}
			}
		#endif
	}
};

#endif