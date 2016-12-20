#ifndef flexDiagonalOperator_H
#define flexDiagonalOperator_H

#include "vector"
#include "flexLinearOperator.h"

template < typename T, typename Tvector >
class flexDiagonalOperator : public flexLinearOperator<T, Tvector>
{
private:
	Tvector diagonalElements;
	T* diagonalElementsPtr;
public:

	flexDiagonalOperator(std::vector<T> _diagonalElements) : flexLinearOperator<T, Tvector>((int)_diagonalElements.size(), (int)_diagonalElements.size(), diagonalOp)
	{
		this->diagonalElements.resize((int)_diagonalElements.size());

		#if __CUDACC__
			thrust::copy(_diagonalElements.begin(), _diagonalElements.end(), this->diagonalElements.begin());
			diagonalElementsPtr = thrust::raw_pointer_cast(this->diagonalElements.data());

		#else
			this->diagonalElements = _diagonalElements;
			this->diagonalElementsPtr = this->diagonalElements.data();
		#endif
	};

	#if __CUDACC__
	flexDiagonalOperator(Tvector _diagonalElements) : diagonalElements(_diagonalElements), flexLinearOperator<T, Tvector>((int)_diagonalElements.size(), (int)_diagonalElements.size(), diagonalOp)
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

	T timesElement(int index, const T* input)
	{
		return this->diagonalElementsPtr[index] * input[index];
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

	#if __CUDACC__
	__device__ T timesElementCUDA(int index, const T* input)
	{
		return this->diagonalElementsPtr[index] * input[index];
	}

	__device__ T getRowsumElementCUDA(int index)
	{
		return fabs(this->diagonalElementsPtr[index]);
	}
	#endif
};

#endif
