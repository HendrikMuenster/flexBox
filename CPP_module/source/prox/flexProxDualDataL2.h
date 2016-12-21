#ifndef flexProxDualL2_H
#define flexProxDualL2_H

#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualDataL2 : public flexProx<T, Tvector>
{
private:

public:

	flexProxDualDataL2() : flexProx<T, Tvector>(dualL2DataProx)
	{
	}

	~flexProxDualDataL2()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		
	}
    
    #if __CUDACC__
        template < typename T >
        struct flexProxDualDataL2Functor
        {
            __host__ __device__
            flexProxDualDataL2Functor(T _alpha) : alpha(_alpha){};

            template <typename Tuple>
            __host__ __device__
            void operator()(Tuple t)
            {
                thrust::get<0>(t) = alpha / (thrust::get<2>(t) + alpha) * (thrust::get<1>(t) - thrust::get<2>(t) * thrust::get<3>(t));
            }

            const T alpha;
        };
    #endif
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{
		#if __CUDACC__
            for (int i = 0; i < dualNumbers.size(); i++)
			{
                auto startIterator = thrust::make_zip_iterator(thrust::make_tuple(data->y[dualNumbers[i]].begin(), data->yTilde[dualNumbers[i]].begin(), data->sigmaElt[dualNumbers[i]].begin(),  fList[i].begin()));
                auto endIterator = thrust::make_zip_iterator(  thrust::make_tuple(data->y[dualNumbers[i]].end(),   data->yTilde[dualNumbers[i]].end(),   data->sigmaElt[dualNumbers[i]].end(),    fList[i].end()));
                
                thrust::for_each(startIterator,endIterator,flexProxDualDataL2Functor<T>(alpha));
            }
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
					ptrY[j] = alpha / (ptrSigma[j] + alpha) * (ptrYtilde[j] - ptrSigma[j] * ptrF[j]);
				}
			}
		#endif
	}
};

#endif