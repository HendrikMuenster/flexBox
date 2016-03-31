#ifndef flexGradientOperator_H
#define flexGradientOperator_H

#include "flexVector.h"
#include "flexLinearOperator.h"

template < typename T >
class flexGradientOperator : public flexLinearOperator<T>
{
private:
	flexVector<int> inputDimension;
	int gradDirection;
	int type;
	bool transposed;

	flexVector<int> gradient2Dxf_1;
	flexVector<int> gradient2Dxf_2;
	flexVector<int> gradient2Dxf_3;

	flexVector<int> gradient2Dxm_1;
	flexVector<int> gradient2Dxm_2;
	flexVector<int> gradient2Dxm_3;

	flexVector<int> gradient2DxfT_1;
	flexVector<int> gradient2DxfT_2;
	flexVector<int> gradient2DxfT_3;

	flexVector<int> gradient2DxmT_1;
	flexVector<int> gradient2DxmT_2;
	flexVector<int> gradient2DxmT_3;

	flexVector<int> gradient2Dyf_1;
	flexVector<int> gradient2Dyf_2;
	flexVector<int> gradient2Dyf_3;

	flexVector<int> gradient2Dym_1;
	flexVector<int> gradient2Dym_2;
	flexVector<int> gradient2Dym_3;

	flexVector<int> gradient2DyfT_1;
	flexVector<int> gradient2DyfT_2;
	flexVector<int> gradient2DyfT_3;

	flexVector<int> gradient2DymT_1;
	flexVector<int> gradient2DymT_2;
	flexVector<int> gradient2DymT_3;
public:

	//type:
	// 0 = forward
	// 1 = backward
	flexGradientOperator(flexVector<int> _inputDimension, int _gradDirection, int _type) : flexLinearOperator<T>(_inputDimension.product(), _inputDimension.product())
	{
		this->inputDimension = _inputDimension;
		this->gradDirection = _gradDirection;
		this->type = _type;
		this->transposed = false;

		//precalculate indices for xForward and yForward gradient
		for (int j = 0; j < this->inputDimension[1]; ++j)
		{
			for (int i = 0; i < this->inputDimension[0]; ++i)
			{
				//x-derivative
				if (i < this->inputDimension[0] - 1)
				{
					// u_x(i,j) forward is u(i+1,j) - u(i,j)
					gradient2Dxf_1.push_back(this->index2DtoLinear(i, j));
					gradient2Dxf_2.push_back(this->index2DtoLinear(i + 1, j));
					gradient2Dxf_3.push_back(this->index2DtoLinear(i, j));
				}
				if (i > 0)
				{
					// u_x(i,j) backward is u(i,j) - u(i-1,j)
					gradient2Dxm_1.push_back(this->index2DtoLinear(i, j));
					gradient2Dxm_2.push_back(this->index2DtoLinear(i, j));
					gradient2Dxm_3.push_back(this->index2DtoLinear(i - 1, j));
				}
				if (i > 0 && i < this->inputDimension[0] - 1)
				{
					// u_x(i,j) forward transposed is u(i-1,j) - u(i,j)
					gradient2DxfT_1.push_back(this->index2DtoLinear(i, j));
					gradient2DxfT_2.push_back(this->index2DtoLinear(i - 1, j));
					gradient2DxfT_3.push_back(this->index2DtoLinear(i, j));
				}
				if (i > 0 && i < this->inputDimension[0] - 1)
				{
					// u_x(i,j) backward transposed is u(i,j) - u(i+1,j)
					gradient2DxmT_1.push_back(this->index2DtoLinear(i, j));
					gradient2DxmT_2.push_back(this->index2DtoLinear(i, j));
					gradient2DxmT_3.push_back(this->index2DtoLinear(i+1, j));
				}

				//y-derivative
				if (j < this->inputDimension[1] - 1)
				{
					// u_y(i,j) forward is u(i,j+1) - u(i,j)
					gradient2Dyf_1.push_back(this->index2DtoLinear(i, j));
					gradient2Dyf_2.push_back(this->index2DtoLinear(i, j + 1));
					gradient2Dyf_3.push_back(this->index2DtoLinear(i, j));
				}
				if (j > 0)
				{
					// u_y(i,j) backward is u(i,j) - u(i,j-1)
					gradient2Dym_1.push_back(this->index2DtoLinear(i, j));
					gradient2Dym_2.push_back(this->index2DtoLinear(i, j));
					gradient2Dym_3.push_back(this->index2DtoLinear(i, j - 1));
				}
				if (j > 0 && j < this->inputDimension[1] - 1)
				{
					// u_y(i,j) forward transposed is u(i,j-1) - u(i,j)
					gradient2DyfT_1.push_back(this->index2DtoLinear(i, j));
					gradient2DyfT_2.push_back(this->index2DtoLinear(i, j - 1));
					gradient2DyfT_3.push_back(this->index2DtoLinear(i, j));
				}
				if (j > 0 && j < this->inputDimension[1] - 1)
				{
					// u_y(i,j) backward transposed is u(i,j) - u(i,j+1)
					gradient2DymT_1.push_back(this->index2DtoLinear(i, j));
					gradient2DymT_2.push_back(this->index2DtoLinear(i, j));
					gradient2DymT_3.push_back(this->index2DtoLinear(i, j + 1));
				}
			}
		}
	};

	flexGradientOperator<T>* copy()
	{
		flexGradientOperator<T>* A = new flexGradientOperator(this->inputDimension, this->gradDirection, this->type);

		return A;
	}

	//apply linear operator to vector
	void times(const flexVector<T> &input, flexVector<T> &output)
	{

	}

	void calculateDifference(const flexVector<T> &input, flexVector<T> &output, flexVector<int> vec1, const flexVector<int> vec2, const flexVector<int> vec3)
	{
		int numElements = vec1.size();
		#pragma omp parallel for
		for (int i = 0; i < numElements; ++i)
		{
			output[vec1[i]] += input[vec2[i]] - input[vec3[i]];
		}
	}

	void timesPlus(const flexVector<T> &input, flexVector<T> &output)
	{
		if (this->inputDimension.size() == 2)
		{
			if (this->gradDirection == 0)
			{
				if (this->type == 0)
				{
					if (this->transposed == false)
					{
						calculateDifference(input, output, gradient2Dxf_1, gradient2Dxf_2, gradient2Dxf_3);
					}
					else
					{
						calculateDifference(input, output, gradient2DxfT_1, gradient2DxfT_2, gradient2DxfT_3);

						//update boundary values differently
						#pragma omp parallel for
						for (int j = 0; j < this->inputDimension[1]; ++j)
						{
							output[this->index2DtoLinear(0, j)] += -input[this->index2DtoLinear(0, j)];
							output[this->index2DtoLinear(this->inputDimension[0] - 1, j)] += input[this->index2DtoLinear(this->inputDimension[0] - 2, j)];
						}
					}
				}
				else if (this->type == 1)
				{
					if (this->transposed == false)
					{
						calculateDifference(input, output, gradient2Dxm_1, gradient2Dxm_2, gradient2Dxm_3);
					}
					else
					{
						calculateDifference(input, output, gradient2DxmT_1, gradient2DxmT_2, gradient2DxmT_3);

						//update boundary values differently
						#pragma omp parallel for
						for (int j = 0; j < this->inputDimension[1]; ++j)
						{
							output[this->index2DtoLinear(0, j)] += -input[this->index2DtoLinear(1, j)];
							output[this->index2DtoLinear(this->inputDimension[0] - 1, j)] += input[this->index2DtoLinear(this->inputDimension[0] - 1, j)];
						}
					}
				}
			}
			else if (this->gradDirection == 1)
			{
				if (this->type == 0)
				{
					if (this->transposed == false)
					{
						calculateDifference(input, output, gradient2Dyf_1, gradient2Dyf_2, gradient2Dyf_3);
					}
					else
					{
						calculateDifference(input, output, gradient2DyfT_1, gradient2DyfT_2, gradient2DyfT_3);

						//update boundary values differently
						#pragma omp parallel for
						for (int i = 0; i < this->inputDimension[0]; ++i)
						{
							output[this->index2DtoLinear(i, 0)] += -input[this->index2DtoLinear(i, 0)];
							output[this->index2DtoLinear(i, this->inputDimension[1] - 1)] += input[this->index2DtoLinear(i, this->inputDimension[1] - 2)];
						}
					}
				}
				else if (this->type == 1)
				{
					if (this->transposed == false)
					{
						calculateDifference(input, output, gradient2Dym_1, gradient2Dym_2, gradient2Dym_3);
					}
					else
					{
						calculateDifference(input, output, gradient2DymT_1, gradient2DymT_2, gradient2DymT_3);

						//update boundary values differently
						#pragma omp parallel for
						for (int i = 0; i < this->inputDimension[0]; ++i)
						{
							output[this->index2DtoLinear(i, 0)] += -input[this->index2DtoLinear(i, 1)];
							output[this->index2DtoLinear(i, this->inputDimension[1] - 1)] += input[this->index2DtoLinear(i, this->inputDimension[1] - 1)];
						}
					}
				}
			}
		}
		else if (this->inputDimension.size() == 3)
		{
			//todo
		}
	}

	void timesMinus(const flexVector<T> &input, flexVector<T> &output)
	{
		if (this->inputDimension.size() == 2)
		{
			if (this->gradDirection == 0)
			{
				if (this->type == 0)
				{
					if (this->transposed == false)
					{
						calculateDifference(input, output, gradient2Dxf_1, gradient2Dxf_3, gradient2Dxf_2);
					}
					else
					{
						calculateDifference(input, output, gradient2DxfT_1, gradient2DxfT_3, gradient2DxfT_2);

						//update boundary values differently
						#pragma omp parallel for
						for (int j = 0; j < this->inputDimension[1]; ++j)
						{
							output[this->index2DtoLinear(0, j)] -= -input[this->index2DtoLinear(0, j)];
							output[this->index2DtoLinear(this->inputDimension[0] - 1, j)] -= input[this->index2DtoLinear(this->inputDimension[0] - 2, j)];
						}
					}
				}
				else if (this->type == 1)
				{
					if (this->transposed == false)
					{
						calculateDifference(input, output, gradient2Dxm_1, gradient2Dxm_3, gradient2Dxm_2);
					}
					else
					{
						calculateDifference(input, output, gradient2DxmT_1, gradient2DxmT_3, gradient2DxmT_2);

						//update boundary values differently
						#pragma omp parallel for
						for (int j = 0; j < this->inputDimension[1]; ++j)
						{
							output[this->index2DtoLinear(0, j)] -= -input[this->index2DtoLinear(1, j)];
							output[this->index2DtoLinear(this->inputDimension[0] - 1, j)] -= input[this->index2DtoLinear(this->inputDimension[0] - 1, j)];
						}
					}
				}
			}
			else if (this->gradDirection == 1)
			{
				if (this->type == 0)
				{
					if (this->transposed == false)
					{
						calculateDifference(input, output, gradient2Dyf_1, gradient2Dyf_3, gradient2Dyf_2);
					}
					else
					{
						calculateDifference(input, output, gradient2DyfT_1, gradient2DyfT_3, gradient2DyfT_2);

						//update boundary values differently
						#pragma omp parallel for
						for (int i = 0; i < this->inputDimension[0]; ++i)
						{
							output[this->index2DtoLinear(i, 0)] -= -input[this->index2DtoLinear(i, 0)];
							output[this->index2DtoLinear(i, this->inputDimension[1] - 1)] -= input[this->index2DtoLinear(i, this->inputDimension[1] - 2)];
						}
					}
				}
				else if (this->type == 1)
				{
					if (this->transposed == false)
					{
						calculateDifference(input, output, gradient2Dym_1, gradient2Dym_3, gradient2Dym_2);
					}
					else
					{
						calculateDifference(input, output, gradient2DymT_1, gradient2DymT_3, gradient2DymT_2);

						//update boundary values differently
						#pragma omp parallel for
						for (int i = 0; i < this->inputDimension[0]; ++i)
						{
							output[this->index2DtoLinear(i, 0)] -= -input[this->index2DtoLinear(i, 1)];
							output[this->index2DtoLinear(i, this->inputDimension[1] - 1)] -= input[this->index2DtoLinear(i, this->inputDimension[1] - 1)];
						}
					}
				}
			}
		}
		else if (this->inputDimension.size() == 3)
		{
			//todo
		}
	}

	T getMaxRowSumAbs()
	{
		//row sum of absolute values is at maximum 2
		return static_cast<T>(2);
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
		return (int)(i + j*this->inputDimension[0]);
	}

	int index3DtoLinear(int i, int j, int k)
	{
		return (int)(i + j*this->inputDimension[0] + k*this->inputDimension[0] * this->inputDimension[1]);
	}
};

#endif
