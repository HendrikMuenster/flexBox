#ifndef flexSuperpixelOperator_H
#define flexSuperpixelOperator_H

#include "vector"
#include "flexLinearOperator.h"

template < typename T, typename Tvector >
class flexSuperpixelOperator : public flexLinearOperator<T, Tvector>
{
private:
	std::vector<int> targetDimension;
	T upsamplingFactor;
	bool transposed;
public:

	flexSuperpixelOperator(std::vector<int> targetDimension_, T upsamplingFactor_) : flexLinearOperator<T, Tvector>((int)(vectorProduct(targetDimension_)), (int)(vectorProduct(targetDimension_)*upsamplingFactor_*upsamplingFactor_), superpixelOp)
	{
		this->targetDimension.resize(targetDimension_.size());
		this->transposed = false;

        this->targetDimension = targetDimension_;
		this->upsamplingFactor = upsamplingFactor_;
	};

	flexSuperpixelOperator<T, Tvector>* copy()
	{
		return new flexSuperpixelOperator<T, Tvector>(this->targetDimension, this->upsamplingFactor);
	}

	int indexI(int index, int sizeX)
	{
		return index % sizeX;
	}

	int indexJ(int index, int sizeX, int sizeY)
	{
		return (index / sizeX) % sizeY;
	}

	int index2DtoLinear(int i, int j, int sizeY)
	{
		return (i*sizeY + j);
	}

	void calcTimes(const Tvector &input, Tvector &output, mySign signRule)
	{

		T factor = (T)1 / (this->upsamplingFactor*this->upsamplingFactor);

		int iOuterEnd = targetDimension[0];
		int jOuterEnd = targetDimension[1];

		int sizeY = targetDimension[1] * (int)this->upsamplingFactor;

		#pragma omp parallel for
		for (int i = 0; i < iOuterEnd; ++i)
		{
			for (int j = 0; j < jOuterEnd; ++j)
			{
				//printf("Output: (%d,%d) : %d\n", i, j, index2DtoLinear(i, j, this->targetDimension[1]));

				int outputIndex = index2DtoLinear(i, j, targetDimension[1]);

				int iInnerStart = i*(int)this->upsamplingFactor;
				int iInnerEnd = (i + 1)*(int)this->upsamplingFactor;

				int jInnerStart = j*(int)this->upsamplingFactor;
				int jInnerEnd = (j + 1)*(int)this->upsamplingFactor;

				T tmpResult = (T)0;

				for (int iInner = iInnerStart; iInner < iInnerEnd; ++iInner)
				{
					for (int jInner = jInnerStart; jInner < jInnerEnd; ++jInner)
					{
						int inputIndex = index2DtoLinear(iInner, jInner, sizeY);

						tmpResult += input[inputIndex];
						/*printf("Inner: (%d,%d) : %d\n", iInner, jInner, inputIndex);

						int innerJ = indexI(inputIndex, this->targetDimension[0] * this->upsamplingFactor);
						int innerI = indexJ(inputIndex, this->targetDimension[0] * this->upsamplingFactor, this->targetDimension[1] * this->upsamplingFactor);

						printf("Back: (%d,%d) \n", innerI, innerJ);

						int backI = innerI / this->upsamplingFactor;
						int backJ = innerJ / this->upsamplingFactor;

						printf("BackInner: (%d,%d) \n", backI, backJ);

						if (backI != i || backJ != j)
						{
							mexErrMsgTxt("PROBLEM!!!\n");
						}*/
					}
				}

				switch (signRule)
				{
					case PLUS:
					{
						output[outputIndex] += factor*tmpResult;
						break;
					}
					case MINUS:
					{
						output[outputIndex] -= factor*tmpResult;
						break;
					}
					case EQUALS:
					{
						output[outputIndex] = factor*tmpResult;
						break;
					}
				}
			}
		}

		//mexErrMsgTxt("Stop!\n");
	}

	void calcTimesTransposed(const Tvector &input, Tvector &output, mySign signRule)
	{
		T factor = (T)1 / (this->upsamplingFactor*this->upsamplingFactor);

		int sizeX = targetDimension[0] * (int)this->upsamplingFactor;
		int sizeY = targetDimension[1] * (int)this->upsamplingFactor;

		#pragma omp parallel for
		for (int i = 0; i < sizeX; ++i)
		{
			for (int j = 0; j < sizeY; ++j)
			{
				int inputIndex = index2DtoLinear(i, j, sizeY);

				//int innerJ = indexI(inputIndex, this->targetDimension[0] * this->upsamplingFactor);
				//int innerI = indexJ(inputIndex, this->targetDimension[0] * this->upsamplingFactor, this->targetDimension[1] * this->upsamplingFactor);

				//printf("Back: (%d,%d) \n", innerI, innerJ);

				int backI = i / (int)this->upsamplingFactor;
				int backJ = j / (int)this->upsamplingFactor;

				
				int outputIndex = index2DtoLinear(backI, backJ, targetDimension[1]);

				//printf("Back: (%d,%d) %d,%d \n", backI, backJ, inputIndex, outputIndex);

				switch (signRule)
				{
					case PLUS:
					{
						output[inputIndex] += factor*input[outputIndex];
						break;
					}
					case MINUS:
					{
						output[inputIndex] -= factor*input[outputIndex];
						break;
					}
					case EQUALS:
					{
						output[inputIndex] = factor*input[outputIndex];
						break;
					}
				}
			}
		}

		//mexErrMsgTxt("Stop!\n");
	}

	//apply linear operator to vector
	void times(const Tvector &input, Tvector &output)
	{
		if (this->transposed)
		{
			calcTimesTransposed(input, output, EQUALS);
		}
		else
		{
			calcTimes(input, output, EQUALS);
		}
	}

	void timesPlus(const Tvector &input, Tvector &output)
	{
		if (this->transposed)
		{
			calcTimesTransposed(input, output, PLUS);
		}
		else
		{
			calcTimes(input, output, PLUS);
		}
	}

	void timesMinus(const Tvector &input, Tvector &output)
	{
		if (this->transposed)
		{
			calcTimesTransposed(input, output, MINUS);
		}
		else
		{
			calcTimes(input, output, MINUS);
		}
	}
    
	std::vector<T> getAbsRowSum()
	{
		if (this->transposed)
		{
			return std::vector<T>(this->getNumRows(), (T)1 / (T)(this->upsamplingFactor*this->upsamplingFactor));
		}
		else
		{
			return std::vector<T>(this->getNumRows(), (T)1);
		}
	}

	T getMaxRowSumAbs()
	{
		if (this->transposed)
		{
			return (T)1 / (T)(this->upsamplingFactor*this->upsamplingFactor);
		}
		else
		{
			return (T)1;
		}
	}

	//transposing the identity does nothing
	void transpose()
	{
		int numRowsTmp = this->getNumRows();
		this->setNumRows(this->getNumCols());
		this->setNumCols(numRowsTmp);

		this->transposed = true;
	}
    
    #ifdef __CUDACC__
    thrust::device_vector<T> getAbsRowSumCUDA()
	{
		Tvector result(this->getNumRows(),(T)1);

		return result;
	}
    #endif
};

#endif
