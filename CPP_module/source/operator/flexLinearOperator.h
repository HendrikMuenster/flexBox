#ifndef flexLinearOperator_H
#define flexLinearOperator_H

#include "vector"
#include "tools.h"

template < typename T, typename Tvector >
class flexLinearOperator
{
private:
	int numRows;
	int numCols;
public:
	linOp type;

	flexLinearOperator(int _numRows, int _numCols)
	{
		numRows = _numRows; 
		numCols = _numCols;
		type = linearOp;
	}

	virtual ~flexLinearOperator()
	{
		if (VERBOSE > 0) printf("Linear operator destructor");
	}

	flexLinearOperator(int _numRows, int _numCols, linOp _type) : type(_type)
	{
		numRows = _numRows;
		numCols = _numCols;
	}

	int getNumCols() const
	{
		return numCols;
	}

	int getNumRows() const
	{
		return numRows;
	}

	void setNumCols(int _numCols)
	{
		numCols = _numCols;
	}

	void setNumRows(int _numRows)
	{
		numRows = _numRows;
	}

	virtual flexLinearOperator<T, Tvector>* copy() = 0;

	//apply linear operator to vector
	virtual void times(const Tvector &input, Tvector &output) = 0;

	virtual void timesPlus(const Tvector &input, Tvector &output) = 0;

	virtual void timesMinus(const Tvector &input, Tvector &output) = 0;

	virtual std::vector<T> getAbsRowSum() = 0;

	#if __CUDACC__		
		virtual thrust::device_vector<T> getAbsRowSumCUDA() = 0;
	#endif

	//used for preconditioning
	virtual T getMaxRowSumAbs() = 0;

	//transpose current matrix
	virtual void transpose() = 0;
};

#endif
