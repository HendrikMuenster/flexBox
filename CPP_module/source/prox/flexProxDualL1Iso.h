#ifndef flexProxDualL1Iso_H
#define flexProxDualL1Iso_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualL1Iso : public flexProx<T, Tvector>
{
private:

public:

    flexProxDualL1Iso() : flexProx<T, Tvector>(dualL1IsoProx){}

	~flexProxDualL1Iso()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}
    
    #ifdef __CUDACC__
	struct flexProxDualL1IsoDim2Functor
	{
		__host__ __device__
		flexProxDualL1IsoDim2Functor(T alpha) : alpha(alpha){}

		template <typename Tuple>
		__host__ __device__
        void operator()(Tuple t)
		{
			T norm = max((T)1, sqrt( pow(thrust::get<2>(t),(int)2) + pow(thrust::get<3>(t),(int)2)) / this->alpha);

			thrust::get<0>(t) = thrust::get<2>(t) / norm;
			thrust::get<1>(t) = thrust::get<3>(t) / norm;
		}

		T alpha;
	};
    #endif

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		#ifdef __CUDACC__
            if (dualNumbers.size() == 2)
			{
                auto startIterator = thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[0]].begin(), data->y[dualNumbers[1]].begin(), data->yTilde[dualNumbers[0]].begin(), data->yTilde[dualNumbers[1]].begin()));
                auto endIterator =   thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[0]].end(),   data->y[dualNumbers[1]].end(),   data->yTilde[dualNumbers[0]].end(),   data->yTilde[dualNumbers[1]].end()));
                
                thrust::for_each(startIterator,endIterator,flexProxDualL1IsoDim2Functor(alpha));
			}
            else
            {
                printf("Alert! Iso prox not implemented in CUDA for dim!=2");
            }
		#else
			if (dualNumbers.size() == 1)
			{
				T* ptrY0 = data->y[dualNumbers[0]].data();

				T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();

				int numElements = (int)data->yTilde[dualNumbers[0]].size();

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					T yTmp = (T)1 / myMax<T>((T)1, fabs(ptrYtilde0[i]) / alpha);

					ptrY0[i] = ptrYtilde0[i] * yTmp;
				}
			}
			else if (dualNumbers.size() == 2)
			{
				T* ptrY0 = data->y[dualNumbers[0]].data();
				T* ptrY1 = data->y[dualNumbers[1]].data();

				T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
				T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();

				int numElements = (int)data->yTilde[dualNumbers[0]].size();

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i]) + pow2(ptrYtilde1[i])) / alpha);

					ptrY0[i] = ptrYtilde0[i] * yTmp;
					ptrY1[i] = ptrYtilde1[i] * yTmp;
				}
			}
			else if (dualNumbers.size() == 3)
			{
				T* ptrY0 = data->y[dualNumbers[0]].data();
				T* ptrY1 = data->y[dualNumbers[1]].data();
				T* ptrY2 = data->y[dualNumbers[2]].data();

				T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
				T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();
				T* ptrYtilde2 = data->yTilde[dualNumbers[2]].data();

				int numElements = (int)data->yTilde[dualNumbers[0]].size();

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i]) + pow2(ptrYtilde1[i]) + pow2(ptrYtilde2[i])) / alpha);

					ptrY0[i] = ptrYtilde0[i] * yTmp;
					ptrY1[i] = ptrYtilde1[i] * yTmp;
					ptrY2[i] = ptrYtilde2[i] * yTmp;
				}
			}
			else
			{
				printf("Alert! Iso prox not implemented for dim>3");
			}
		#endif
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{

	}
};

#endif