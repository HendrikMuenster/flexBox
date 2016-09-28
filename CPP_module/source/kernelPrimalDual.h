


#ifndef kernelPrimalDual_H
#define kernelPrimalDual_H

#include "tools.h"
#include "flexBox.h"
#include "flexTermDual.h"
#include "flexBoxDataGPU.h"
#include "flexDiagonalOperator.h"
#include "flexLinearOperator.h"
#include "flexIdentityOperator.h"
#include "flexGradientOperator.h"
#include "flexMatrixGPU.h"
#include "flexSuperpixelOperator.h"

#include "flexDualizedDataTerm.h"

#include <thrust/device_vector.h> 

// Could be any number, but the whole array should fit into shared memory 

template<typename T, typename Tvector>
__global__ void updateY(int numTerms, int numDualVars, flexTermDual<T, Tvector>** termDualList, T** ptrYList, T** ptrYOldList, T** prtXBar, T** sigmaList, int** dualNumbers, int** primalNumbers, const int* numElementsY)
{
	const int index = threadIdx.x + blockIdx.x * blockDim.x;

	T yTildeList[CONST_ARRAY_SIZE];
	for (int dualNum = 0; dualNum < numDualVars; ++dualNum)
	{
		if (index < numElementsY[dualNum])
		{
			yTildeList[dualNum] = (T)0;
		}
	}
	
	for (int termNumber = 0; termNumber < numTerms; ++termNumber)
	{
		int* dualNumbersTerm = dualNumbers[termNumber];
		int* primalNumbersTerm = primalNumbers[termNumber];

		auto* termDual = termDualList[termNumber];

		int numDualsTerm = termDual->getNumberDuals();//  numDuals[termNumber];
		int numPrimalsTerm = termDual->getNumberPrimals();//numPrimals[termNumber];

		for (int i = 0; i < numDualsTerm; ++i)
		{
			const int dualNum = dualNumbersTerm[i];

			if (index < numElementsY[dualNum])
			{
				for (int j = 0; j < numPrimalsTerm; ++j)
				{
					int operatorNumber = numPrimalsTerm * i + j;

					const int primalNum = primalNumbersTerm[j];

					switch (termDual->operatorListG[operatorNumber]->type)
					{
					case diagonalOp:
						{
							yTildeList[dualNum] += static_cast<flexDiagonalOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXBar[primalNum]);
							break;
						}
						case identityOp:
						{
							yTildeList[dualNum] += static_cast<flexIdentityOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXBar[primalNum]);
							break;
						}
						case gradientOp:
						{
							yTildeList[dualNum] += static_cast<flexGradientOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXBar[primalNum]);
							break;
						}
						case matrixGPUOp:
						{
							yTildeList[dualNum] += static_cast<flexMatrixGPU<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXBar[primalNum]);
							break;
						}
						case superpixelOp:
						{
							yTildeList[dualNum] += static_cast<flexSuperpixelOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXBar[primalNum]);
							break;
						}
					}
				}
			}
		}
	}

	for (int dualNum = 0; dualNum < numDualVars; ++dualNum)
	{
		if (index < numElementsY[dualNum])
		{
			yTildeList[dualNum] = ptrYOldList[dualNum][index] + sigmaList[dualNum][index] * yTildeList[dualNum];
		}
	}


	for (int termNumber = 0; termNumber < numTerms; ++termNumber)
	{
		int* dualNumbersTerm = dualNumbers[termNumber];

		//toDo: change this and shift into applyProxElement
		if (index < numElementsY[dualNumbersTerm[0]])
		{
			int* primalNumbersTerm = primalNumbers[termNumber];

			flexTermDual<T, Tvector>* termDual = termDualList[termNumber];

			int numDualsTerm = termDual->getNumberDuals();//  numDuals[termNumber];
			int numPrimalsTerm = termDual->getNumberPrimals();//numPrimals[termNumber];
			//do prox
			switch (termDual->p)
			{
				case dualL1AnisoProx:
				case dualL1IsoProx:
				case dualL2Prox:
				case dualFrobeniusProx:
				case dualHuberProx:
				{
					static_cast<flexDualizedOperator<T, Tvector>*>(termDual)->applyProxElement(ptrYList, yTildeList, sigmaList, dualNumbersTerm, numDualsTerm, index);
					break;
				}
				case dualL2DataProx:
				case dualL1DataProx:
				case dualKLDataProx:
				{
					static_cast<flexDualizedDataTerm<T, Tvector>*>(termDual)->applyProxElement(ptrYList, yTildeList, sigmaList, dualNumbersTerm, numDualsTerm, index);
					break;
				}
			}
		}
	}
};


template<typename T, typename Tvector>
__global__ void errorY(int numTerms, int numDualVars, flexTermDual<T, Tvector>** termDualList, T** ptrYList, T** ptrYOldList, T** ptrYError, T** prtX, T** prtXOld, T** sigmaList, int** dualNumbers, int** primalNumbers,const int* numElementsY)
{
	const int index = threadIdx.x + blockIdx.x * blockDim.x;

	for (int dualNum = 0; dualNum < numDualVars; ++dualNum)
	{
		if (index < numElementsY[dualNum])
		{
			ptrYError[dualNum][index] = fabs(ptrYList[dualNum][index] - ptrYOldList[dualNum][index]) / sigmaList[dualNum][index];
		}
	}

	//go through all dual terms
	for (int termNumber = 0; termNumber < numTerms; ++termNumber)
	{
		int* dualNumbersTerm = dualNumbers[termNumber];
		int* primalNumbersTerm = primalNumbers[termNumber];

		flexTermDual<T, Tvector>* termDual = termDualList[termNumber];

		int numDualsTerm = termDual->getNumberDuals();//  numDuals[termNumber];
		int numPrimalsTerm = termDual->getNumberPrimals();//numPrimals[termNumber];

		for (int i = 0; i < numDualsTerm; ++i)
		{
			const int dualNum = dualNumbersTerm[i];

			if (index < numElementsY[dualNum])
			{
				//load yTilde from memory
				T yErr = (T)0;

				for (int j = 0; j < numPrimalsTerm; ++j)
				{
					int operatorNumber = numPrimalsTerm * i + j;

					const int primalNum = primalNumbersTerm[j];

					switch (termDual->operatorListG[operatorNumber]->type)
					{
						case diagonalOp:
						{
							yErr += static_cast<flexDiagonalOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtX[primalNum])
								- static_cast<flexDiagonalOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXOld[primalNum]);
							break;
						}
						case identityOp:
						{
							yErr += static_cast<flexIdentityOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtX[primalNum])
								- static_cast<flexIdentityOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXOld[primalNum]);
							break;
						}
						case gradientOp:
						{
							yErr += static_cast<flexGradientOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtX[primalNum])
								- static_cast<flexGradientOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXOld[primalNum]);
							break;
						}
						case matrixGPUOp:
						{
							yErr += static_cast<flexMatrixGPU<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtX[primalNum])
								- static_cast<flexMatrixGPU<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXOld[primalNum]);
							break;
						}
						case superpixelOp:
						{
							yErr += static_cast<flexSuperpixelOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtX[primalNum])
								- static_cast<flexSuperpixelOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->timesElementCUDA(index, prtXOld[primalNum]);
							break;
						}
					}
				}
				//save yTilde to local memory
				ptrYError[dualNum][index] -= yErr;
			}
		}
	}
};


template<typename T, typename Tvector>
__global__ void errorX(int numTermsDual, int numDualVars, int numPrimalVars, flexTermDual<T, Tvector>** termDualList, T** ptrYList, T** ptrYOldList, T** ptrXError, T** prtX, T** prtXOld, T** tauList, int** dualNumbers, int** primalNumbers, int* numElementsX)
{
	const int index = threadIdx.x + blockIdx.x * blockDim.x;

	for (int primalNum = 0; primalNum < numPrimalVars; ++primalNum)
	{
		if (index < numElementsX[primalNum])
		{
			ptrXError[primalNum][index] = fabs(prtX[primalNum][index] - prtXOld[primalNum][index]) / tauList[primalNum][index];
		}
	}

	for (int termNumber = 0; termNumber < numTermsDual; ++termNumber)
	{
		int* dualNumbersTerm = dualNumbers[termNumber];
		int* primalNumbersTerm = primalNumbers[termNumber];

		flexTermDual<T, Tvector>* termDual = termDualList[termNumber];

		int numDualsTerm = termDual->getNumberDuals();
		int numPrimalsTerm = termDual->getNumberPrimals();

		for (int i = 0; i < numDualsTerm; ++i)
		{
			const int dualNum = dualNumbersTerm[i];

			for (int j = 0; j < numPrimalsTerm; ++j)
			{
				int operatorNumber = numPrimalsTerm * i + j;

				int primalNum = primalNumbersTerm[j];

				if (index < numElementsX[primalNum])
				{
					switch (termDual->operatorListTG[operatorNumber]->type)
					{
						case diagonalOp:
						{
							ptrXError[primalNum][index] += -static_cast<flexDiagonalOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYOldList[dualNum]) +
								static_cast<flexDiagonalOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
						case identityOp:
						{
							ptrXError[primalNum][index] += -static_cast<flexIdentityOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYOldList[dualNum]) +
								static_cast<flexIdentityOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
						case gradientOp:
						{
							ptrXError[primalNum][index] += -static_cast<flexGradientOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYOldList[dualNum]) +
								static_cast<flexGradientOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
						case matrixGPUOp:
						{
							ptrXError[primalNum][index] += -static_cast<flexMatrixGPU<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYOldList[dualNum]) +
								static_cast<flexMatrixGPU<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
						case superpixelOp:
						{
							ptrXError[primalNum][index] += -static_cast<flexSuperpixelOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYOldList[dualNum]) +
								static_cast<flexSuperpixelOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
					}
				}
			}
		}
	}



};

//template<typename T, typename Tvector, int numPrimals>
template<typename T, typename Tvector>
__global__ void updateX(int numTermsDual, int numTermsPrimal, int numPrimalVars, flexTermDual<T, Tvector>** termDualList, flexTermPrimal<T, Tvector>** termPrimalList, T** ptrXList, T** ptrXOldList, T** ptrXTildeList, T** ptrXBarList, T** ptrYList, T** tauList, int** dualNumbers, int** primalNumbers, int** primalCorrespPrimal, int* numElementsX)
{
	//x value corresponds to index in vector
	//y value corresponds to variable number

	const int index = threadIdx.x + blockIdx.x * blockDim.x;

	__shared__ int s_numElementsX[CONST_ARRAY_SIZE];

	if (threadIdx.x == 0)
	{
		for (int primalNum = 0; primalNum < numPrimalVars; ++primalNum)
		{
			s_numElementsX[primalNum] = numElementsX[primalNum];
		}
	}

	syncthreads();

	T xTildeList[CONST_ARRAY_SIZE];

	for (int primalNum = 0; primalNum < numPrimalVars; ++primalNum)
	{
		if (index < s_numElementsX[primalNum])
		{
			xTildeList[primalNum] = (T)0;
		}
	}

	
	//calculate xTilde using all dual terms
	for (int termNumber = 0; termNumber < numTermsDual; ++termNumber)
	{
		int* dualNumbersTerm = dualNumbers[termNumber];
		int* primalNumbersTerm = primalNumbers[termNumber];

		flexTermDual<T, Tvector>* termDual = termDualList[termNumber];

		int numDualsTerm = termDual->getNumberDuals();
		int numPrimalsTerm = termDual->getNumberPrimals();

		for (int i = 0; i < numDualsTerm; ++i)
		{
			const int dualNum = dualNumbersTerm[i];

			for (int j = 0; j < numPrimalsTerm; ++j)
			{
				int operatorNumber = numPrimalsTerm * i + j;

				int primalNum = primalNumbersTerm[j];

				if (index < s_numElementsX[primalNum])
				{
					switch (termDual->operatorListTG[operatorNumber]->type)
					{
						case diagonalOp:
						{
							xTildeList[primalNum] += static_cast<flexDiagonalOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
						case identityOp:
						{
							xTildeList[primalNum] += static_cast<flexIdentityOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
						case gradientOp:
						{
							xTildeList[primalNum] += static_cast<flexGradientOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
						case matrixGPUOp:
						{
							xTildeList[primalNum] += static_cast<flexMatrixGPU<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
						case superpixelOp:
						{
							xTildeList[primalNum] += static_cast<flexSuperpixelOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->timesElementCUDA(index, ptrYList[dualNum]);
							break;
						}
					}

					
				}
			}
		}
	}

	for (int primalNum = 0; primalNum < numPrimalVars; ++primalNum)
	{
		if (index < s_numElementsX[primalNum])
		{
			xTildeList[primalNum] = ptrXOldList[primalNum][index] - tauList[primalNum][index] * xTildeList[primalNum];
		}
	}

	//primal prox
	for (int termNumber = 0; termNumber < numTermsPrimal; ++termNumber)
	{
		flexTermPrimal<T, Tvector>* termPrimal = termPrimalList[termNumber];

		int numPrimalsTerm = termPrimal->getNumberPrimals();

		if (index < s_numElementsX[primalCorrespPrimal[termNumber][0]])
		{
			switch (termPrimal->p)
			{
				case primalEmptyProx:
				{
					static_cast<flexTermPrimal<T, Tvector>*>(termPrimal)->applyProxElement(ptrXList, xTildeList, tauList, primalCorrespPrimal[termNumber], numPrimalsTerm, index);
					break;
				}
			}
		}
	}

	//overrelaxation
	for (int primalNum = 0; primalNum < numPrimalVars; ++primalNum)
	{
		if (index < s_numElementsX[primalNum])
		{
			ptrXBarList[primalNum][index] = 2 * xTildeList[primalNum] - ptrXOldList[primalNum][index];
		}
	}
};

typedef float T;
typedef thrust::device_vector<T> Tvector;

template __global__ void updateX<T, Tvector>(int, int, int, flexTermDual<T, Tvector>**, flexTermPrimal<T, Tvector>**, T**, T**, T**, T**, T**, T**, int**, int**, int**, int*);
//template __global__ void updateX<T, Tvector, 2>(int, int, int, flexTermDual<T, Tvector>**, flexTermPrimal<T, Tvector>**, T**, T**, T**, T**, T**, T**, int**, int**, int**, int*);
//template __global__ void updateX<T, Tvector, 3>(int, int, int, flexTermDual<T, Tvector>**, flexTermPrimal<T, Tvector>**, T**, T**, T**, T**, T**, T**, int**, int**, int**, int*);

template<typename T, typename Tvector>
__global__ void sigmaElement(int numTerms, int numDualVars, flexTermDual<T, Tvector>** termDualList, T** sigmaList, int** dualNumbers, int** primalNumbers, int* numElementsY)
{
	const int index = threadIdx.x + blockIdx.x * blockDim.x;

	for (int dualNum = 0; dualNum < numDualVars; ++dualNum)
	{
		if (index < numElementsY[dualNum])
		{
			sigmaList[dualNum][index] = (T)0;
		}
	}

	for (int termNumber = 0; termNumber < numTerms; ++termNumber)
	{
		int* dualNumbersTerm = dualNumbers[termNumber];
		int* primalNumbersTerm = primalNumbers[termNumber];

		flexTermDual<T, Tvector>* termDual = termDualList[termNumber];

		int numDualsTerm = termDual->getNumberDuals();//  numDuals[termNumber];
		int numPrimalsTerm = termDual->getNumberPrimals();//numPrimals[termNumber];

		for (int i = 0; i < numDualsTerm; ++i)
		{
			const int dualNum = dualNumbersTerm[i];

			if (index < numElementsY[dualNum])
			{
				//load yTilde from memory
				for (int j = 0; j < numPrimalsTerm; ++j)
				{
					int operatorNumber = numPrimalsTerm * i + j;

					const int primalNum = primalNumbersTerm[j];
					switch (termDual->operatorListG[operatorNumber]->type)
					{
						case diagonalOp:
						{
							sigmaList[dualNum][index] += static_cast<flexDiagonalOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
						case identityOp:
						{
							sigmaList[dualNum][index] += static_cast<flexIdentityOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
						case gradientOp:
						{
							sigmaList[dualNum][index] += static_cast<flexGradientOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
						case matrixGPUOp:
						{
							sigmaList[dualNum][index] += static_cast<flexMatrixGPU<T, Tvector>*>(termDual->operatorListG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
						case superpixelOp:
						{
							sigmaList[dualNum][index] += static_cast<flexSuperpixelOperator<T, Tvector>*>(termDual->operatorListG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
					}
				}
			}
		}
	}

	for (int dualNum = 0; dualNum < numDualVars; ++dualNum)
	{
		if (index < numElementsY[dualNum])
		{
			sigmaList[dualNum][index] = (T)1 / myMaxGPUf((T)0.001, sigmaList[dualNum][index]);
		}
	}
}

template<typename T, typename Tvector>
__global__ void tauElement(int numTerms, int numPrimalVars, flexTermDual<T, Tvector>** termDualList, T** tauList, int** dualNumbers, int** primalNumbers, int* numElementsX)
{
	const int index = threadIdx.x + blockIdx.x * blockDim.x;

	for (int primalNum = 0; primalNum < numPrimalVars; ++primalNum)
	{
		if (index < numElementsX[primalNum])
		{
			tauList[primalNum][index] = (T)0;
		}
	}

	for (int termNumber = 0; termNumber < numTerms; ++termNumber)
	{
		int* dualNumbersTerm = dualNumbers[termNumber];
		int* primalNumbersTerm = primalNumbers[termNumber];

		flexTermDual<T, Tvector>* termDual = termDualList[termNumber];

		int numDualsTerm = termDual->getNumberDuals();
		int numPrimalsTerm = termDual->getNumberPrimals();

		for (int i = 0; i < numDualsTerm; ++i)
		{
			const int dualNum = dualNumbersTerm[i];

			for (int j = 0; j < numPrimalsTerm; ++j)
			{
				int operatorNumber = numPrimalsTerm * i + j;

				int primalNum = primalNumbersTerm[j];

				if (index < numElementsX[primalNum])
				{
					switch (termDual->operatorListTG[operatorNumber]->type)
					{
						case diagonalOp:
						{
							tauList[primalNum][index] += static_cast<flexDiagonalOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
						case identityOp:
						{
							tauList[primalNum][index] += static_cast<flexIdentityOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
						case gradientOp:
						{
							tauList[primalNum][index] += static_cast<flexGradientOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
						case matrixGPUOp:
						{
							tauList[primalNum][index] += static_cast<flexMatrixGPU<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
						case superpixelOp:
						{
							tauList[primalNum][index] += static_cast<flexSuperpixelOperator<T, Tvector>*>(termDual->operatorListTG[operatorNumber])->getRowsumElementCUDA(index);
							break;
						}
					}
				}
			}
		}
	}

	for (int primalNum = 0; primalNum < numPrimalVars; ++primalNum)
	{
		if (index < numElementsX[primalNum])
		{
			tauList[primalNum][index] = (T)1 / myMaxGPUf((T)0.001, tauList[primalNum][index]);
		}
	}
	
};


#endif