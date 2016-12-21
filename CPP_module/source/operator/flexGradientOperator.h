#ifndef flexGradientOperator_H
#define flexGradientOperator_H

#include "vector"
#include "flexLinearOperator.h"

template < typename T, typename Tvector >
class flexGradientOperator : public flexLinearOperator<T, Tvector>
{
private:
	std::vector<int> inputDimension;
	
	int* inputDimensionPtr;
	int gradDirection;
	int type;
	bool transposed;
	int numberDimensions;

public:

	//type:
	// 0 = forward
	// 1 = backward
	flexGradientOperator(std::vector<int> _inputDimension, int _gradDirection, int _type) : flexLinearOperator<T, Tvector>(vectorProduct(_inputDimension), vectorProduct(_inputDimension),gradientOp)
	{
		
		this->gradDirection = _gradDirection;
		this->type = _type;
		this->transposed = false;
		this->numberDimensions = (int)_inputDimension.size();

		this->inputDimension = _inputDimension;
		this->inputDimensionPtr = this->inputDimension.data();

	};

	flexGradientOperator<T, Tvector>* copy()
	{

		std::vector<int> dimsCopy;
		dimsCopy.resize(this->inputDimension.size());

		std::copy(this->inputDimension.begin(), this->inputDimension.end(), dimsCopy.begin());

		return new flexGradientOperator<T, Tvector>(dimsCopy, this->gradDirection, this->type);
	}

	//apply linear operator to vector
	void times(const Tvector &input, Tvector &output)
	{

	}



	void dxp2d(const Tvector &input, Tvector &output, mySign s)
	{
		int sizeY = this->inputDimension[1];
		int sizeX = this->inputDimension[0] - 1;

		#pragma omp parallel for
		for (int j = 0; j < sizeY; ++j)
		{
			for (int i = 0; i < sizeX; ++i)
			{
				const int tmpIndex = this->index2DtoLinear(i, j);

				switch (s)
				{
					case PLUS:
					{
						output[tmpIndex] += input[tmpIndex + 1] - input[tmpIndex];
						break;
					}
					case MINUS:
					{
						output[tmpIndex] -= input[tmpIndex + 1] - input[tmpIndex];
						break;
					}
					case EQUALS:
					{
						output[tmpIndex] = input[tmpIndex + 1] - input[tmpIndex];
						break;
					}
				}
			}
		}
	}

	void dyp2d(const Tvector &input, Tvector &output, mySign s)
	{
		int sizeY = this->inputDimension[1] - 1;
		int sizeX = this->inputDimension[0];

		#pragma omp parallel for
		for (int j = 0; j < sizeY; ++j)
		{
			for (int i = 0; i < sizeX; ++i)
			{
				const int tmpIndex = this->index2DtoLinear(i, j);

				switch (s)
				{
					case PLUS:
					{
						output[tmpIndex] += input[tmpIndex + sizeX] - input[tmpIndex];
						break;
					}
					case MINUS:
					{
						output[tmpIndex] -= input[tmpIndex + sizeX] - input[tmpIndex];
						break;
					}
					case EQUALS:
					{
						output[tmpIndex] = input[tmpIndex + sizeX] - input[tmpIndex];
						break;
					}
				}
			}
		}
	}

	void dxp2dTransposed(const Tvector &input, Tvector &output, mySign s)
	{
		int sizeY = this->inputDimension[1];
		int sizeX = this->inputDimension[0] - 1;

		#pragma omp parallel for
		for (int j = 0; j < sizeY; ++j)
		{
			for (int i = 1; i < sizeX; ++i)
			{
				int tmpIndex = this->index2DtoLinear(i, j);

				switch (s)
				{
					case PLUS:
					{
						output[tmpIndex] += -(input[tmpIndex] - input[tmpIndex - 1]);
						break;
					}
					case MINUS:
					{
						output[tmpIndex] -= -(input[tmpIndex] - input[tmpIndex - 1]);
						break;
					}
					case EQUALS:
					{
						output[tmpIndex] = -(input[tmpIndex] - input[tmpIndex - 1]);
						break;
					}
				}
			}
		}

		for (int j = 0; j < this->inputDimension[1]; ++j)
		{
			switch (s)
			{
				case PLUS:
				{
					output[this->index2DtoLinear(0, j)] += -input[this->index2DtoLinear(0, j)];
					output[this->index2DtoLinear(this->inputDimension[0] - 1, j)] += input[this->index2DtoLinear(this->inputDimension[0] - 2, j)];
					break;
				}
				case MINUS:
				{
					output[this->index2DtoLinear(0, j)] -= -input[this->index2DtoLinear(0, j)];
					output[this->index2DtoLinear(this->inputDimension[0] - 1, j)] -= input[this->index2DtoLinear(this->inputDimension[0] - 2, j)];
					break;
				}
				case EQUALS:
				{
					output[this->index2DtoLinear(0, j)] = -input[this->index2DtoLinear(0, j)];
					output[this->index2DtoLinear(this->inputDimension[0] - 1, j)] = input[this->index2DtoLinear(this->inputDimension[0] - 2, j)];
					break;
				}
			}
		}
	}

	void dyp2dTransposed(const Tvector &input, Tvector &output, mySign s)
	{
		int sizeY = this->inputDimension[1] - 1;
		int sizeX = this->inputDimension[0];

		#pragma omp parallel for
		for (int j = 1; j < sizeY; ++j)
		{
			for (int i = 0; i < sizeX; ++i)
			{
				int tmpIndex = this->index2DtoLinear(i, j);

				switch (s)
				{
					case PLUS:
					{
						output[tmpIndex] += -(input[tmpIndex] - input[tmpIndex - sizeX]);
						break;
					}
					case MINUS:
					{
						output[tmpIndex] -= -(input[tmpIndex] - input[tmpIndex - sizeX]);
						break;
					}
					case EQUALS:
					{
						output[tmpIndex] = -(input[tmpIndex] - input[tmpIndex - sizeX]);
						break;
					}
				}
			}
		}

		for (int i = 0; i < this->inputDimension[0]; ++i)
		{
			switch (s)
			{
			case PLUS:
			{
				output[this->index2DtoLinear(i, 0)] += -input[this->index2DtoLinear(i, 0)];
				output[this->index2DtoLinear(i, this->inputDimension[1] - 1)] += input[this->index2DtoLinear(i, this->inputDimension[1] - 2)];
				break;
			}
			case MINUS:
			{
				output[this->index2DtoLinear(i, 0)] -= -input[this->index2DtoLinear(i, 0)];
				output[this->index2DtoLinear(i, this->inputDimension[1] - 1)] -= input[this->index2DtoLinear(i, this->inputDimension[1] - 2)];
				break;
			}
			case EQUALS:
			{
				output[this->index2DtoLinear(i, 0)] = -input[this->index2DtoLinear(i, 0)];
				output[this->index2DtoLinear(i, this->inputDimension[1] - 1)] = input[this->index2DtoLinear(i, this->inputDimension[1] - 2)];
				break;
			}
			}
		}
	}


	void timesPlus(const Tvector &input, Tvector &output)
	{
        #if __CUDACC__
			dim3 block2d = dim3(32, 16, 1);
			dim3 grid2d = dim3((this->inputDimension[0] + block2d.x - 1) / block2d.x, (this->inputDimension[1] + block2d.y - 1) / block2d.y, 1);

			T* ptrOutput = thrust::raw_pointer_cast(output.data());
			const T* ptrInput = thrust::raw_pointer_cast(input.data());
		#endif

		if (this->inputDimension.size() == 2)
		{
			if (this->gradDirection == 0)
			{
				if (this->transposed == false)
				{
					#if __CUDACC__
						dxp2dCUDA << <grid2d, block2d >> >(ptrOutput,ptrInput, this->inputDimension[0], this->inputDimension[1], SIGN_PLUS);
					#else
						this->dxp2d(input,output,PLUS);
					#endif
				}
				else
				{
					#if __CUDACC__
						dxp2dTransposedCUDA << <grid2d, block2d >> >(ptrOutput, ptrInput, this->inputDimension[0], this->inputDimension[1], SIGN_PLUS);
					#else
						this->dxp2dTransposed(input, output, PLUS);
					#endif
				}
			}
			else if (this->gradDirection == 1)
			{
				if (this->transposed == false)
				{
					#if __CUDACC__
						dyp2dCUDA << <grid2d, block2d >> >(ptrOutput, ptrInput, this->inputDimension[0], this->inputDimension[1], SIGN_PLUS);
					#else
						this->dyp2d(input, output, PLUS);
					#endif
				}
				else
				{
					#if __CUDACC__
						dyp2dTransposedCUDA << <grid2d, block2d >> >(ptrOutput, ptrInput, this->inputDimension[0], this->inputDimension[1], SIGN_PLUS);
					#else
						this->dyp2dTransposed(input, output, PLUS);
					#endif
				}
			}
		}
		else
		{
            printf("Gradient not implemented for dim!=2");
			//todo
		}
	}

	void timesMinus(const Tvector &input, Tvector &output)
	{
		#if __CUDACC__
			dim3 block2d = dim3(32, 16, 1);
			dim3 grid2d = dim3((this->inputDimension[0] + block2d.x - 1) / block2d.x, (this->inputDimension[1] + block2d.y - 1) / block2d.y, 1);

			T* ptrOutput = thrust::raw_pointer_cast(output.data());
			const T* ptrInput = thrust::raw_pointer_cast(input.data());
		#endif
        
		if (this->inputDimension.size() == 2)
		{
			if (this->gradDirection == 0)
			{
				if (this->transposed == false)
				{
					#if __CUDACC__
						dxp2dCUDA << <grid2d, block2d >> >(ptrOutput, ptrInput, this->inputDimension[0], this->inputDimension[1], SIGN_MINUS);
					#else
						this->dxp2d(input, output, MINUS);
					#endif
				}
				else
				{
					#if __CUDACC__
						dxp2dTransposedCUDA << <grid2d, block2d >> >(ptrOutput, ptrInput, this->inputDimension[0], this->inputDimension[1], SIGN_MINUS);
					#else
						this->dxp2dTransposed(input, output, MINUS);
					#endif
				}
			}
			else if (this->gradDirection == 1)
			{
				if (this->transposed == false)
				{
					#if __CUDACC__
						dyp2dCUDA << <grid2d, block2d >> >(ptrOutput, ptrInput, this->inputDimension[0], this->inputDimension[1], SIGN_MINUS);
					#else
						this->dyp2d(input, output, MINUS);
					#endif
				}
				else
				{
					#if __CUDACC__
						dyp2dTransposedCUDA << <grid2d, block2d >> >(ptrOutput, ptrInput, this->inputDimension[0], this->inputDimension[1], SIGN_MINUS);
					#else
						this->dyp2dTransposed(input, output, MINUS);
					#endif
				}
			}
		}
		else
        {
            printf("Gradient not implemented for dim!=2");
			//todo
		}
	}

	T getMaxRowSumAbs() const
	{
		//row sum of absolute values is at maximum 2
		return static_cast<T>(2);
	}

	T getMaxRowSumAbs()
	{
		//row sum of absolute values is at maximum 2
		return static_cast<T>(2);
	}

	std::vector<T> getAbsRowSum()
	{
		std::vector<T> result(this->getNumRows(),(T)2);

		return result;
	}

	//transposing the identity does nothing
	void transpose()
	{
		//swith type from forward to backward and vice versa
		if (this->transposed == false)
		{
			this->transposed = true;
		}
		else
		{
			this->transposed = false;
		}
	}

	int index2DtoLinear(int i, int j)
	{
		return (i + j*this->inputDimension[0]);
	}

	#if __CUDACC__	
	thrust::device_vector<T> getAbsRowSumCUDA()
	{
		thrust::device_vector<T> result(this->getNumRows(), (T)2);

		return result;
	}
	#endif
};

#endif
