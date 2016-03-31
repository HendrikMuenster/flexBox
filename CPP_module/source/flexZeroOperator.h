#ifndef flexZeroOperator_H
#define flexZeroOperator_H

#include "flexVector.h"
#include "flexLinearOperator.h"

template < typename T >
class flexZeroOperator : public flexLinearOperator<T>
{
public:
	flexZeroOperator(int _numRows, int _numCols) : flexLinearOperator<T>(_numRows, _numCols){};

	flexZeroOperator<T>* copy()
	{
		flexZeroOperator<T>* A = new flexZeroOperator(this->getNumRows(), this->getNumCols());

		return A;
	}


	//apply linear operator to vector
	void times(const flexVector<T> &input, flexVector<T> &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < this->getNumRows(); ++i)
		{
			output[i] = static_cast<T>(0);
		}
	}

	void timesPlus(const flexVector<T> &input, flexVector<T> &output)
	{
	}

	void timesMinus(const flexVector<T> &input, flexVector<T> &output)
	{
	}

	T getMaxRowSumAbs()
	{
		return static_cast<T>(1);
	}

	//transposing the identity does nothing
	void transpose()
	{
		int numRowsTmp = this->getNumRows();
		this->setNumRows(this->getNumCols());
		this->setNumCols(numRowsTmp);
	}
};

#endif
