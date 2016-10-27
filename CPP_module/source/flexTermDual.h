#ifndef flexTermDual_H
#define flexTermDual_H

//#include "flexProxList.h"
#include "flexBoxData.h"
#include "vector"
#include "flexLinearOperator.h"

#if __CUDACC__
	#include "flexMatrixGPU.h"
#endif

//template < typename T > class flexBoxData;

template < typename T, typename Tvector >
class flexTermDual
{
	private:
		int numberVars;
		int numberPrimals;
		float paramFloat1;

	public:
		const prox p;
		T alpha;
		std::vector<flexLinearOperator<T, Tvector>* > operatorList;
		std::vector<flexLinearOperator<T, Tvector>* > operatorListT;

		#if __CUDACC__
			//list of pointers to operators on GPU
			thrust::device_vector<flexLinearOperator<T, Tvector>* > thrust_operatorList;
			thrust::device_vector<flexLinearOperator<T, Tvector>* > thrust_operatorListT;

			flexLinearOperator<T, Tvector>** operatorListG;
			flexLinearOperator<T, Tvector>** operatorListTG;
		#endif

			flexTermDual(prox _p, T _alpha, int _numberPrimals, int _numberVars, float _paramFloat1) : alpha(_alpha), numberPrimals(_numberPrimals), numberVars(_numberVars), p(_p), paramFloat1(_paramFloat1){}

		virtual ~flexTermDual()
		{ 
			if (VERBOSE > 0) printf("Destructor virtual\n!");
		
			for (int i = operatorList.size() - 1; i >= 0; --i)
			{
#if __CUDACC__
				//free operator memory
				cudaFree(thrust_operatorList[i]); CUDA_CHECK;
				cudaFree(thrust_operatorListT[i]); CUDA_CHECK;
#endif
				delete operatorList[i];
				delete operatorListT[i];
			}

			operatorList.clear();
			operatorListT.clear();
		};


		int getNumberVars()
		{
			return numberVars;
		}

		int dualVarLength(int num)
		{
			return this->operatorList[num]->getNumRows();
		}
		
		virtual void applyProx(flexBoxData<T, Tvector> *data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers) = 0;

#if __CUDACC__
		__device__ int getNumberPrimals()
		{
			return this->numberPrimals;
		}

		__device__ int getNumberDuals()
		{
			return this->numberVars;
		}

		__device__ void applyProxElement(T yUpdate[CONST_ARRAY_SIZE], const T yTilde[CONST_ARRAY_SIZE], T** sigmaList, const int* dualNumbers, int numDuals, int index)
		{
			return (T)0;
		}

		void createGPUOperators()
		{
			//create operators on GPU

			thrust_operatorList.resize(this->operatorList.size());
			thrust_operatorListT.resize(this->operatorList.size());

			for (int i = 0; i < this->operatorList.size(); ++i)
			{
				flexLinearOperator<T, Tvector>* operatorPointer;
				flexLinearOperator<T, Tvector>* operatorPointerT;

				//copy operators to GPU
				switch (this->operatorList[i]->type)
				{
					case diagonalOp:
					{
						cudaMalloc((void **)&operatorPointer, sizeof(flexDiagonalOperator<T, Tvector>)); CUDA_CHECK;
						cudaMalloc((void **)&operatorPointerT, sizeof(flexDiagonalOperator<T, Tvector>)); CUDA_CHECK;
						cudaMemcpy(operatorPointer, this->operatorList[i], sizeof(flexDiagonalOperator<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						cudaMemcpy(operatorPointerT, this->operatorListT[i], sizeof(flexDiagonalOperator<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						break;
					}
					case identityOp:
					{
						cudaMalloc((void **)&operatorPointer, sizeof(flexIdentityOperator<T, Tvector>)); CUDA_CHECK;
						cudaMalloc((void **)&operatorPointerT, sizeof(flexIdentityOperator<T, Tvector>)); CUDA_CHECK;
						cudaMemcpy(operatorPointer, this->operatorList[i], sizeof(flexIdentityOperator<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						cudaMemcpy(operatorPointerT, this->operatorListT[i], sizeof(flexIdentityOperator<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						break;
					}
					case gradientOp:
					{
						cudaMalloc((void **)&operatorPointer, sizeof(flexGradientOperator<T, Tvector>)); CUDA_CHECK;
						cudaMalloc((void **)&operatorPointerT, sizeof(flexGradientOperator<T, Tvector>)); CUDA_CHECK;
						cudaMemcpy(operatorPointer, this->operatorList[i], sizeof(flexGradientOperator<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						cudaMemcpy(operatorPointerT, this->operatorListT[i], sizeof(flexGradientOperator<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						break;
					}
					case matrixGPUOp:
					{
						cudaMalloc((void **)&operatorPointer, sizeof(flexMatrixGPU<T, Tvector>)); CUDA_CHECK;
						cudaMalloc((void **)&operatorPointerT, sizeof(flexMatrixGPU<T, Tvector>)); CUDA_CHECK;
						cudaMemcpy(operatorPointer, this->operatorList[i], sizeof(flexMatrixGPU<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						cudaMemcpy(operatorPointerT, this->operatorListT[i], sizeof(flexMatrixGPU<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						break;
					}
					case superpixelOp:
					{
						cudaMalloc((void **)&operatorPointer, sizeof(flexSuperpixelOperator<T, Tvector>)); CUDA_CHECK;
						cudaMalloc((void **)&operatorPointerT, sizeof(flexSuperpixelOperator<T, Tvector>)); CUDA_CHECK;
						cudaMemcpy(operatorPointer, this->operatorList[i], sizeof(flexSuperpixelOperator<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						cudaMemcpy(operatorPointerT, this->operatorListT[i], sizeof(flexSuperpixelOperator<T, Tvector>), cudaMemcpyHostToDevice); CUDA_CHECK;
						break;
					}
				}

				thrust_operatorList[i] = operatorPointer;
				thrust_operatorListT[i] = operatorPointerT;
			}

			operatorListG = thrust::raw_pointer_cast(thrust_operatorList.data());
			operatorListTG = thrust::raw_pointer_cast(thrust_operatorListT.data());
		}
#endif
};

#endif