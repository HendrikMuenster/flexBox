#ifndef flexLinearOperator_H
#define flexLinearOperator_H

#include "flexVector.h"

template < typename T >
class flexLinearOperator
{
private:
	int numRows;
	int numCols;
public:
	flexLinearOperator(int _numRows, int _numCols)
	{
		numRows = _numRows; 
		numCols = _numCols;
	}

	int getNumCols()
	{
		return numCols;
	}

	int getNumRows()
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

	virtual flexLinearOperator<T>* copy() = 0;

	//apply linear operator to vector
	virtual void times(const flexVector<T> &input, flexVector<T> &output) = 0;

	virtual void timesPlus(const flexVector<T> &input, flexVector<T> &output) = 0;

	virtual void timesMinus(const flexVector<T> &input, flexVector<T> &output) = 0;

	//used for preconditioning
	virtual T getMaxRowSumAbs() = 0;

	//transpose current matrix
	virtual void transpose() = 0;
};

#endif
