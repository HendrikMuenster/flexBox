#ifndef flexProxDualHuber_H
#define flexProxDualHuber_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualHuber : public flexProx<T, Tvector>
{
private:
	T huberEpsilon;
public:

	flexProxDualHuber(T _huberEpsilon) : flexProx<T, Tvector>(dualHuberProx)
	{
		huberEpsilon = _huberEpsilon;
	}

	~flexProxDualHuber()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}
    
    #if __CUDACC__
	struct flexProxDualHuberDim2Functor
	{
		__host__ __device__
		flexProxDualHuberDim2Functor(T _epsi, T _alpha) : epsi(_epsi), alpha(_alpha), epsiAlpha(epsi / alpha){}

		template <typename Tuple>
		__host__ __device__
        void operator()(Tuple t)
		{
            T huberFactor1 = (T)1 / ((T)1 + thrust::get<4>(t) * epsiAlpha);
            T huberFactor2 = (T)1 / ((T)1 + thrust::get<5>(t) * epsiAlpha);
            
			T norm = max((T)1, sqrt( pow(thrust::get<2>(t)*huberFactor1,(int)2) + pow(thrust::get<3>(t)*huberFactor2,(int)2)) / alpha);

			thrust::get<0>(t) = thrust::get<2>(t) * huberFactor1 / norm;
			thrust::get<1>(t) = thrust::get<3>(t) * huberFactor2 / norm;
		}
        
        const T epsi;
		const T alpha;
        const T epsiAlpha;
	};
    #endif

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		#if __CUDACC__
            if (dualNumbers.size() == 2)
			{
                auto startIterator = thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[0]].begin(), data->y[dualNumbers[1]].begin(), data->yTilde[dualNumbers[0]].begin(), data->yTilde[dualNumbers[1]].begin(), data->sigmaElt[dualNumbers[0]].begin(), data->sigmaElt[dualNumbers[1]].begin()));
                auto endIterator =   thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[0]].end(),   data->y[dualNumbers[1]].end(),   data->yTilde[dualNumbers[0]].end(),   data->yTilde[dualNumbers[1]].end(),   data->sigmaElt[dualNumbers[0]].end(),   data->sigmaElt[dualNumbers[1]].end()));
                
                thrust::for_each(startIterator,endIterator,flexProxDualHuberDim2Functor(this->huberEpsilon,alpha));
			}
            else
            {
                printf("Alert! Huber prox not implemented in CUDA for dim!=2\n");
            }
		#else
			if (dualNumbers.size() == 1)
			{
				T* ptrY0 = data->y[dualNumbers[0]].data();

				T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();

				T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

				int numElements = (int)data->yTilde[dualNumbers[0]].size();
                
                T epsiAlpha = this->huberEpsilon / alpha;

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * epsiAlpha);

					T yTmp = (T)1 / myMax<T>((T)1, fabs(ptrYtilde0[i] * huberFactor) / alpha);

					ptrY0[i] = ptrYtilde0[i] * huberFactor * yTmp;
				}
			}
			else if (dualNumbers.size() == 2)
			{
				T* ptrY0 = data->y[dualNumbers[0]].data();
				T* ptrY1 = data->y[dualNumbers[1]].data();

				T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
				T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();

				T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

				int numElements = (int)data->yTilde[dualNumbers[0]].size();
                
                T epsiAlpha = this->huberEpsilon / alpha;

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * epsiAlpha);

					T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i] * huberFactor) + pow2(ptrYtilde1[i] * huberFactor)) / alpha);

					ptrY0[i] = ptrYtilde0[i] * huberFactor * yTmp;
					ptrY1[i] = ptrYtilde1[i] * huberFactor * yTmp;
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

				T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

				int numElements = (int)data->yTilde[dualNumbers[0]].size();
                
                T epsiAlpha = this->huberEpsilon / alpha;

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * epsiAlpha);

					T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i] * huberFactor) + pow2(ptrYtilde1[i] * huberFactor) + pow2(ptrYtilde2[i] * huberFactor)) / alpha);

					ptrY0[i] = ptrYtilde0[i] * huberFactor * yTmp;
					ptrY1[i] = ptrYtilde1[i] * huberFactor * yTmp;
					ptrY2[i] = ptrYtilde2[i] * huberFactor * yTmp;
				}
			}
			else
			{
				printf("Alert! Huber prox not implemented for dim>3");
			}
		#endif
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{

	}
};

#endif