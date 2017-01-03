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
    
    #if __CUDACC__
	struct flexProxDualDataL1Functor
	{
		__host__ __device__
		flexProxDualDataL1Functor(T alpha) : alpha(alpha){}

		template <typename Tuple>
		__host__ __device__
        void operator()(Tuple t)
		{
			thrust::get<0>(t) = min(this->alpha, max(-this->alpha,thrust::get<1>(t) - alpha * thrust::get<2>(t) * thrust::get<3>(t)));
		}

		T alpha;
	};
    #endif

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{
		#if __CUDACC__
            for (int k = 0; k < dualNumbers.size(); k++)
            {
                auto startIterator = thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[k]].begin(), data->yTilde[dualNumbers[k]].begin(), data->sigmaElt[dualNumbers[k]].begin(), fList[k].begin()));
                auto endIterator =   thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[k]].end(),   data->yTilde[dualNumbers[k]].end(),   data->sigmaElt[dualNumbers[k]].end(),   fList[k].end()));
                
                thrust::for_each(startIterator,endIterator,flexProxDualDataL1Functor(alpha));
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
					ptrY[j] = myMin<T>(alpha, myMax<T>(-alpha, ptrYtilde[j] - alpha * ptrSigma[j] * ptrF[j]));
				}
			}
		#endif
	}
};

#endif