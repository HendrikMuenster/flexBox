#ifndef flexProxDualBoxConstraint_H
#define flexProxDualBoxConstraint_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualBoxConstraint : public flexProx<T, Tvector>
{
private:
	T minVal;
	T maxVal;
public:

	flexProxDualBoxConstraint(T _minVal, T _maxVal) : flexProx<T, Tvector>(dualBoxConstraintProx)
	{
		minVal = _minVal;
		maxVal = _maxVal;
	}

	~flexProxDualBoxConstraint()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}
    
    #if __CUDACC__
	struct flexProxDualBoxConstraintFunctor
	{
		__host__ __device__
		flexProxDualBoxConstraintFunctor(T _minVal, T _maxVal) : minVal(_minVal), maxVal(_maxVal){}

		template <typename Tuple>
		__host__ __device__
        void operator()(Tuple t)
		{
			thrust::get<0>(t) = min((T)0, thrust::get<1>(t) - thrust::get<2>(t) * this->minVal) + max((T)0, thrust::get<1>(t) - thrust::get<2>(t) * this->maxVal);
		}

		const T minVal;
        const T maxVal;
	};
    #endif

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		#if __CUDACC__
            for (int k = 0; k < dualNumbers.size(); k++)
            {
                auto startIterator = thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[k]].begin(), data->yTilde[dualNumbers[k]].begin(), data->sigmaElt[dualNumbers[k]].begin()));
                auto endIterator =   thrust::make_zip_iterator( thrust::make_tuple(data->y[dualNumbers[k]].end(),   data->yTilde[dualNumbers[k]].end(),   data->sigmaElt[dualNumbers[k]].end()));
                
                thrust::for_each(startIterator,endIterator,flexProxDualBoxConstraintFunctor(this->minVal,this->maxVal));
            }
		#else
			for (int k = 0; k < dualNumbers.size(); k++)
			{
				T* ptrY = data->y[dualNumbers[k]].data();
				T* ptrYtilde = data->yTilde[dualNumbers[k]].data();
				T* ptrSigma = data->sigmaElt[dualNumbers[k]].data();

				int numElements = (int)data->yTilde[dualNumbers[k]].size();

				#pragma omp parallel for
				for (int i = 0; i < numElements; i++)
				{
					ptrY[i] = myMax<T>((T)0, ptrYtilde[i] - ptrSigma[i] * this->maxVal) + myMin<T>((T)0, ptrYtilde[i] - ptrSigma[i] * this->minVal);
				}
			}
		#endif
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<Tvector> &fList)
	{

	}
};

#endif