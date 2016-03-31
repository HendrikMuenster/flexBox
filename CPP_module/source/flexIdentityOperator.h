#ifndef flexIdentityOperator_H
#define flexIdentityOperator_H

#include "flexVector.h"
#include "flexLinearOperator.h"

template < typename T >
class flexIdentityOperator : public flexLinearOperator<T>
{
private:
	bool minus;
public:

	flexIdentityOperator(int _numRows, int _numCols,bool _minus) : flexLinearOperator<T>(_numRows, _numCols)
	{
		minus = _minus;
	};

	flexIdentityOperator<T>* copy()
	{
		flexIdentityOperator<T>* A = new flexIdentityOperator(this->getNumRows(), this->getNumCols(),this->minus);

		return A;
	}

	//apply linear operator to vector
	void times(const flexVector<T> &input, flexVector<T> &output)
	{
		int numElements = output.size();

		if (this->minus == true)
		{
			for (int i = 0; i < numElements; ++i)
			{
				output[i] -= input[i];
			}
		}
		else
		{
			for (int i = 0; i < numElements; ++i)
			{
				output[i] = input[i];
			}
		}
	}

	void doTimesPlus(const flexVector<T> &input, flexVector<T> &output)
	{
		int numElements = output.size();

		for (int i = 0; i < numElements; ++i)
		{
			output[i] += input[i];
		}
	}

	void doTimesMinus(const flexVector<T> &input, flexVector<T> &output)
	{
		int numElements = output.size();
		for (int i = 0; i < numElements; ++i)
		{
			output[i] -= input[i];
		}
	}

	void timesPlus(const flexVector<T> &input, flexVector<T> &output)
	{
		if (this->minus == true)
		{
			doTimesMinus(input, output);
		}
		else
		{
			doTimesPlus(input, output);
		}
	}

	void timesMinus(const flexVector<T> &input, flexVector<T> &output)
	{
		if (this->minus == true)
		{
			doTimesPlus(input, output);
		}
		else
		{
			doTimesMinus(input, output);
		}
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
