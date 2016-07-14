#ifndef flexDiagonalOperator_H
#define flexDiagonalOperator_H

#include "vector"
#include "flexLinearOperator.h"

template < typename T, typename Tvector >
class flexDiagonalOperator : public flexLinearOperator<T, Tvector>
{
private:
	Tvector diagonalElements;

	#if __CUDACC__
		//save pointer to device memory
		T* diagonalElementsPtr;
	#endif
public:

	flexDiagonalOperator(std::vector<T> _diagonalElements) : flexLinearOperator<T, Tvector>(_diagonalElements.size(), _diagonalElements.size(), diagonalOp)
	{
		this->diagonalElements.resize(_diagonalElements.size());

		#if __CUDACC__
			thrust::copy(_diagonalElements.begin(), _diagonalElements.end(), this->diagonalElements.begin());
			diagonalElementsPtr = thrust::raw_pointer_cast(this->diagonalElements.data());

		#else
			this->diagonalElements = _diagonalElements;
		#endif
	};

	#if __CUDACC__
	flexDiagonalOperator(Tvector _diagonalElements) : diagonalElements(_diagonalElements), flexLinearOperator<T, Tvector>(_diagonalElements.size(), _diagonalElements.size(),diagonalOp)
	{
		diagonalElementsPtr = thrust::raw_pointer_cast(this->diagonalElements.data());
	};
	#endif

	flexDiagonalOperator<T, Tvector>* copy()
	{
		flexDiagonalOperator<T, Tvector>* A = new flexDiagonalOperator<T, Tvector>(this->diagonalElements);

		return A;
	}

	//apply linear operator to vector
	void times(const Tvector &input, Tvector &output)
	{
		vectorAddVectorTimesVector(output, input, this->diagonalElements, SIGN_EQUALS);
	}

	void timesPlus(const Tvector &input, Tvector &output)
	{
		vectorAddVectorTimesVector(output, input, this->diagonalElements, SIGN_PLUS);
	}

	void timesMinus(const Tvector &input, Tvector &output)
	{
		vectorAddVectorTimesVector(output, input, this->diagonalElements, SIGN_MINUS);
	}

	T getMaxRowSumAbs()
	{
		Tvector diagonalElementsCopy = this->diagonalElements;

		vectorAbs(diagonalElementsCopy);

		return vectorMax(diagonalElementsCopy);
	}

	//transposing the identity does nothing
	void transpose(){}

	#if __CUDACC__
	__device__ T timesElement(int index, const T* input)
	{
		return this->diagonalElementsPtr[index] * input[index];
	}

	__device__ T getRowsumElement(int index)
	{
		return fabs(this->diagonalElementsPtr[index]);
	}
	#endif
};

#endif
