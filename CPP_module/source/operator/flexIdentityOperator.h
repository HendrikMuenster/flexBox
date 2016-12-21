#ifndef flexIdentityOperator_H
#define flexIdentityOperator_H


#include "vector"
#include "tools.h"
#include "flexLinearOperator.h"

template < typename T, typename Tvector >
class flexIdentityOperator : public flexLinearOperator<T, Tvector>
{
private:
	bool minus;
public:

	flexIdentityOperator(int _numRows, int _numCols, bool _minus) : flexLinearOperator<T, Tvector>(_numRows, _numCols, identityOp)
	{
		minus = _minus;
	};

	flexIdentityOperator<T, Tvector>* copy()
	{
		flexIdentityOperator<T, Tvector>* A = new flexIdentityOperator<T, Tvector>(this->getNumRows(), this->getNumCols(), this->minus);

		return A;
	}

	//apply linear operator to vector
	void times(const Tvector &input, Tvector &output)
	{
		int numElements = (int)output.size();

		if (this->minus == true)
		{
			int numElements = (int)input.size();
			#pragma omp parallel for
			for (int i = 0; i < numElements; ++i)
			{
				output[i] = -input[i];
			}
		}
		else
		{
			int numElements = (int)input.size();
			#pragma omp parallel for
			for (int i = 0; i < numElements; ++i)
			{
				output[i] = input[i];
			}
		}
	}

	void doTimesPlus(const Tvector &input, Tvector &output)
	{
		#if __CUDACC__
			thrust::transform(output.begin(), output.end(), input.begin(), output.begin(), thrust::plus<T>());
		#else
            int numElements = (int)input.size();
            #pragma omp parallel for
            for (int i = 0; i < numElements; ++i)
            {
                output[i] += input[i];
            }
		#endif
	}

	void doTimesMinus(const Tvector &input, Tvector &output)
	{
		#if __CUDACC__
			thrust::transform(output.begin(), output.end(), input.begin(), output.begin(), thrust::minus<T>());
		#else
            int numElements = (int)input.size();
            #pragma omp parallel for
            for (int i = 0; i < numElements; ++i)
            {
                output[i] -= input[i];
            }
		#endif
	}

	void timesPlus(const Tvector &input, Tvector &output)
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

	void timesMinus(const Tvector &input, Tvector &output)
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

	std::vector<T> getAbsRowSum()
	{
		std::vector<T> result(this->getNumRows(), (T)1);

		return result;
	}

	//transposing the identity does nothing
	void transpose()
	{
		int numRowsTmp = this->getNumRows();
		this->setNumRows(this->getNumCols());
		this->setNumCols(numRowsTmp);
	}

	#if __CUDACC__
	thrust::device_vector<T> getAbsRowSumCUDA()
	{
		thrust::device_vector<T> result(this->getNumRows(), (T)1);

		return result;
	}
	#endif
};

#endif
