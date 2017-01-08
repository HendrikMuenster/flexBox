#ifndef flexProxDualInnerProduct_H
#define flexProxDualInnerProduct_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualInnerProduct : public flexProx<T, Tvector>
{
private:

public:

	flexProxDualInnerProduct() : flexProx<T, Tvector>(dualInnerProductProx)
	{
	}

	~flexProxDualInnerProduct()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}
    


	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{
		#ifdef __CUDACC__
            for (int i = 0; i < dualNumbers.size(); i++)
            {
                thrust::transform(fList[i].begin(), fList[i].end(), data->y[dualNumbers[i]].begin(), alpha * _1);
            }
		#else
			for (int i = 0; i < dualNumbers.size(); i++)
			{
				T* ptrY = data->y[dualNumbers[i]].data();
				T* ptrF = fList[i].data();

				int numElements = (int)data->yTilde[dualNumbers[i]].size();

				#pragma omp parallel for
				for (int j = 0; j < numElements; j++)
				{
					ptrY[j] = alpha * ptrF[j];
				}
			}
		#endif
	}
};

#endif