#ifndef flexDiagonalOperator_H
#define flexDiagonalOperator_H

#include "flexVector.h"
#include "flexLinearOperator.h"

template < typename T >
class flexDiagonalOperator : public flexLinearOperator<T>
{
private:
	flexVector<T> diagonalElements;
public:

	flexDiagonalOperator(flexVector<T> _diagonalElements) : flexLinearOperator(_diagonalElements.size(), _diagonalElements.size())
	{
		this->diagonalElements = _diagonalElements;
	};

	flexDiagonalOperator<T>* copy()
	{
		flexDiagonalOperator<T>* A = new flexDiagonalOperator(this->diagonalElements);

		return A;
	}

	//apply linear operator to vector
	void times(const flexVector<T> &input, flexVector<T> &output)
	{
		int numElements = output.size();

		#pragma omp parallel for
		for (int i = 0; i < numElements; ++i)
		{
			output[i] = this->diagonalElements[i] * input[i];
		}
	}

	void timesPlus(const flexVector<T> &input, flexVector<T> &output)
	{
		int numElements = output.size();

		#pragma omp parallel for
		for (int i = 0; i < numElements; ++i)
		{
			output[i] += this->diagonalElements[i] * input[i];
		}
	}

	void timesMinus(const flexVector<T> &input, flexVector<T> &output)
	{
		int numElements = output.size();

		#pragma omp parallel for
		for (int i = 0; i < numElements; ++i)
		{
			output[i] -= this->diagonalElements[i] * input[i];
		}
	}

	T getMaxRowSumAbs()
	{
		flexVector<T> diagonalElementsCopy = this->diagonalElements;
		diagonalElementsCopy.abs();
		return diagonalElementsCopy.max();
	}

	//transposing the identity does nothing
	void transpose(){}
};

#endif
