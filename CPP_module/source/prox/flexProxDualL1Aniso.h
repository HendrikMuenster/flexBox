#ifndef flexProxDualL1Aniso_H
#define flexProxDualL1Aniso_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualL1Aniso : public flexProx<T, Tvector>
{
private:

public:

    flexProxDualL1Aniso() : flexProx<T, Tvector>(dualL1AnisoProx){}
    
	~flexProxDualL1Aniso()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}
    
    #if __CUDACC__
    template<typename T>
	struct flexProxDualL1AnisoFunctor
	{
		__host__ __device__
		flexProxDualL1AnisoFunctor(T alpha) : alpha(alpha){}

		template <typename Tuple>
		__host__ __device__
        void operator()(Tuple t)
		{
			thrust::get<0>(t) = std::min(this->alpha, std::max(-this->alpha,thrust::get<1>(t)));
		}

		T alpha;
	};
    #endif
    
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		#if __CUDACC__
            for (int k = 0; k < dualNumbers.size(); k++)
            {
                auto startIterator = thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[k]].begin(), data->yTilde[dualNumbers[k]].begin()));
                auto endIterator =   thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[k]].end(),   data->yTilde[dualNumbers[k]].end()));
                
                thrust::for_each(startIterator,endIterator,flexProxDualL1AnisoFunctor<T>(alpha));
            }
		#else
			for (int k = 0; k < dualNumbers.size(); k++)
			{
				T* ptrY = data->y[dualNumbers[k]].data();
				T* ptrYtilde = data->yTilde[dualNumbers[k]].data();

				int numElements = (int)data->yTilde[dualNumbers[k]].size();

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					ptrY[i] = myMin<T>(alpha, myMax<T>(-alpha, ptrYtilde[i]));
				}
			}
		#endif
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{

	}
};

#endif