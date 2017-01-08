#ifndef flexZeroOperator_H
#define flexZeroOperator_H

#include "vector"
#include "flexLinearOperator.h"

template < typename T, typename Tvector >
class flexZeroOperator : public flexLinearOperator<T, Tvector>
{
public:
	flexZeroOperator(int _numRows, int _numCols) : flexLinearOperator<T, Tvector>(_numRows, _numCols, zeroOp){};

	flexZeroOperator<T, Tvector>* copy()
	{
		flexZeroOperator<T, Tvector>* A = new flexZeroOperator<T, Tvector>(this->getNumRows(), this->getNumCols());

		return A;
	}


	//apply linear operator to vector
	void times(const Tvector &input, Tvector &output)
	{
		vectorScalarSet(output, (T)0);
	}

	void timesPlus(const Tvector &input, Tvector &output){}

	void timesMinus(const Tvector &input, Tvector &output){}

	T getMaxRowSumAbs()
	{
		return static_cast<T>(1);
	}

	std::vector<T> getAbsRowSum()
	{
		std::vector<T> result(this->getNumRows(),(T)0);

		return result;
	}

	//transposing the identity does nothing
	void transpose()
	{
		int numRowsTmp = this->getNumRows();
		this->setNumRows(this->getNumCols());
		this->setNumCols(numRowsTmp);
	}

#ifdef __CUDACC__
	thrust::device_vector<T> getAbsRowSumCUDA()
	{
		Tvector result(this->getNumRows(),(T)0);

		return result;
	}
#endif
};

#endif
