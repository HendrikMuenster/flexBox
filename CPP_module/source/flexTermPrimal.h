#ifndef flexTermPrimal_H
#define flexTermPrimal_H

#include "flexProxList.h"
#include "flexBoxData.h"
#include "vector"
#include "tools.h"

template < typename T, typename Tdata>
class flexTermPrimal
{
	private:
		int numberVars;
		
	public:
		T alpha;
		prox p;

		flexTermPrimal(int _numberVars, T _alpha, prox _p) : p(_p)
		{
			this->numberVars = _numberVars;
			this->alpha = _alpha;
		};

		int getNumberVars()
		{
			return numberVars;
		}

		std::vector<int> getDims()
		{
			return this->dims;
		}

		virtual void applyProx(flexBoxData<T, Tdata>* data, std::vector<T> tau, std::vector<int> primalNumbers)
		{
			switch (p)
			{
				case primalEmptyProx :
				{
					flexProxList<T, Tdata>::primalEmptyProx(data, tau, primalNumbers);
					break;
				}
			}
		}

#if __CUDACC__
		__device__ int getNumberPrimals()
		{
			return this->numberVars;
		}

		__device__ void applyProxElement(T** ptrXList,T* xTilde, T** tau, const int* primalNumbers, int numPrimals, int index)
		{
			#pragma unroll
			for (int i = 0; i < numPrimals; ++i)
			{
				const int primalNum = primalNumbers[i];

				ptrXList[primalNum][index] = xTilde[primalNum];
			}
		}
#endif
};

#endif