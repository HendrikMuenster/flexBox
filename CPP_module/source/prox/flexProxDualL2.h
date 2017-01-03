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
    
    #if __CUDACC__
        struct flexProxDualL2Functor
        {
            __host__ __device__
            flexProxDualL2Functor(T _alpha) : alpha(_alpha){};

            template <typename Tuple>
            __host__ __device__
            void operator()(Tuple t)
            {
                thrust::get<0>(t) = this->alpha / (thrust::get<2>(t) + this->alpha) * thrust::get<1>(t);
            }

            const T alpha;
        };
    #endif

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		#if __CUDACC__
            for (int i = 0; i < dualNumbers.size(); i++)
			{
                auto startIterator = thrust::make_zip_iterator(thrust::make_tuple(data->y[dualNumbers[i]].begin(), data->yTilde[dualNumbers[i]].begin(), data->sigmaElt[dualNumbers[i]].begin()));
                auto endIterator = thrust::make_zip_iterator(  thrust::make_tuple(data->y[dualNumbers[i]].end(),   data->yTilde[dualNumbers[i]].end(),   data->sigmaElt[dualNumbers[i]].end()));
                
                thrust::for_each(startIterator,endIterator,flexProxDualL2Functor(alpha));
            }
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
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{

	}
};

#endif