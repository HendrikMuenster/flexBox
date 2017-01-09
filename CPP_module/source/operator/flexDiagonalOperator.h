#ifndef flexDiagonalOperator_H
#define flexDiagonalOperator_H

#include "vector"
#include "flexLinearOperator.h"

template < typename T, typename Tvector >
class flexDiagonalOperator : public flexLinearOperator<T, Tvector>
{
private:
	Tvector diagonalElements;
public:

	flexDiagonalOperator(std::vector<T> _diagonalElements) : flexLinearOperator<T, Tvector>((int)_diagonalElements.size(), (int)_diagonalElements.size(), diagonalOp)
	{
		this->diagonalElements.resize((int)_diagonalElements.size());

		#ifdef __CUDACC__
			thrust::copy(_diagonalElements.begin(), _diagonalElements.end(), this->diagonalElements.begin());

		#else
			this->diagonalElements = _diagonalElements;
		#endif
	};

	#ifdef __CUDACC__
	flexDiagonalOperator(Tvector _diagonalElements) : diagonalElements(_diagonalElements), flexLinearOperator<T, Tvector>((int)_diagonalElements.size(), (int)_diagonalElements.size(), diagonalOp){};
	#endif

	flexDiagonalOperator<T, Tvector>* copy()
	{
		flexDiagonalOperator<T, Tvector>* A = new flexDiagonalOperator<T, Tvector>(this->diagonalElements);

		return A;
	}
    
    #ifdef __CUDACC__
    struct flexDiagonalOperatorTimesFunctor
	{
		__host__ __device__ 
		flexDiagonalOperatorTimesFunctor(){}

		template <typename Tuple>
		__host__ __device__
		void operator()(Tuple t)
		{
            thrust::get<0>(t) = thrust::get<1>(t) * thrust::get<2>(t);
		}
	};
    #endif
    
	//apply linear operator to vector
	void times(const Tvector &input, Tvector &output)
	{
        #ifdef __CUDACC__
            thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(output.begin(), input.begin(), this->diagonalElements.begin())),
                thrust::make_zip_iterator(thrust::make_tuple(output.end(),   input.end(),   this->diagonalElements.end())),
			flexDiagonalOperatorTimesFunctor());
        #else
            int numElements = (int)output.size();

            #pragma omp parallel for
            for (int i = 0; i < numElements; ++i)
            {
                output[i] = input[i] * this->diagonalElements[i];
            }
        #endif
	}
    
    #ifdef __CUDACC__
    struct flexDiagonalOperatorTimesPlusFunctor
	{
		__host__ __device__ 
		flexDiagonalOperatorTimesPlusFunctor(){}

		template <typename Tuple>
		__host__ __device__
		void operator()(Tuple t)
		{
            thrust::get<0>(t) += thrust::get<1>(t) * thrust::get<2>(t);
		}
	};
    #endif
	void timesPlus(const Tvector &input, Tvector &output)
	{
        #ifdef __CUDACC__
            thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(output.begin(), input.begin(), this->diagonalElements.begin())),
                thrust::make_zip_iterator(thrust::make_tuple(output.end(),   input.end(),   this->diagonalElements.end())),
			flexDiagonalOperatorTimesPlusFunctor());
        #else
            int numElements = (int)output.size();

            #pragma omp parallel for
            for (int i = 0; i < numElements; ++i)
            {
                output[i] += input[i] * this->diagonalElements[i];
            }
        #endif
	}

    #ifdef __CUDACC__
    struct flexDiagonalOperatorTimesMinusFunctor
	{
		__host__ __device__ 
		flexDiagonalOperatorTimesMinusFunctor(){}

		template <typename Tuple>
		__host__ __device__
		void operator()(Tuple t)
		{
            thrust::get<0>(t) -= thrust::get<1>(t) * thrust::get<2>(t);
		}
	};
    #endif
	void timesMinus(const Tvector &input, Tvector &output)
	{
        #ifdef __CUDACC__
            thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(output.begin(), input.begin(), this->diagonalElements.begin())),
                thrust::make_zip_iterator(thrust::make_tuple(output.end(),   input.end(),   this->diagonalElements.end())),
			flexDiagonalOperatorTimesPlusFunctor());
        #else
            int numElements = (int)output.size();

            #pragma omp parallel for
            for (int i = 0; i < numElements; ++i)
            {
                output[i] -= input[i] * this->diagonalElements[i];
            }
        #endif
	}

	std::vector<T> getAbsRowSum()
	{
		std::vector<T> result(this->getNumRows());
		#pragma omp parallel for
		for (int k = 0; k < this->getNumRows(); ++k)
		{
			result[k] = std::abs(this->diagonalElements[k]);
		}

		return result;
	}

	T getMaxRowSumAbs()
	{
		Tvector diagonalElementsCopy = this->diagonalElements;

		vectorAbs(diagonalElementsCopy);

		return vectorMax(diagonalElementsCopy);
	}

	//transposing the identity does nothing
	void transpose(){}

	#ifdef __CUDACC__
	thrust::device_vector<T> getAbsRowSumCUDA()
	{
		Tvector diagonalElementsCopy = this->diagonalElements;

		vectorAbs(diagonalElementsCopy);

		return diagonalElementsCopy;
	}
	#endif
};

#endif
