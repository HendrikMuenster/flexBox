#ifndef flexProxDualL1Aniso_H
#define flexProxDualL1Aniso_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualL1Aniso : public flexProx<T, Tvector>
{
private:

public:

	flexProxDualL1Aniso() : flexProx<T, Tvector>()
	{
		
	}

	~flexProxDualL1Aniso()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		for (int k = 0; k < dualNumbers.size(); k++)
		{
			T* ptrY = data->y[dualNumbers[k]].data();
			T* ptrYtilde = data->yTilde[dualNumbers[k]].data();

			int numElements = data->yTilde[dualNumbers[k]].size();

			#pragma omp parallel for
			for (int i = 0; i < numElements; i++)
			{
				ptrY[i] = myMin<T>(alpha, myMax<T>(-alpha, ptrYtilde[i]));
			}
		}
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<std::vector<T>> &fList)
	{

	}
};

#endif