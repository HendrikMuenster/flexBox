#ifndef flexDualizedOperator_H
#define flexDualizedOperator_H

#include "flexTermDual.h"

template < typename T, typename Tvector >
class flexDualizedOperator : public flexTermDual<T, Tvector>
{
public:
	flexDualizedOperator(prox _p, T _alpha, int numberPrimals, std::vector<flexLinearOperator<T, Tvector>* > _operatorList) : flexTermDual<T, Tvector>(_p, _alpha, numberPrimals, _operatorList.size() / numberPrimals)
	{
			
		this->operatorList = _operatorList;

		//create sigma and tau
		this->myTau.resize(numberPrimals, 0.0);
		this->mySigma.resize(_operatorList.size() / numberPrimals, 0.0);

		for (int i = 0; i < _operatorList.size() / numberPrimals; ++i)
		{
			for (int j = 0; j < numberPrimals; ++j)
			{
				int opNum = i*numberPrimals + j;

				this->operatorListT.push_back(this->operatorList[opNum]->copy());
				this->operatorListT[opNum]->transpose();

				this->mySigma[i] += this->operatorList[opNum]->getMaxRowSumAbs();
				this->myTau[j] += this->operatorListT[opNum]->getMaxRowSumAbs();
			}
		}

#if __CUDACC__
		this->createGPUOperators();
#endif
	};

	~flexDualizedOperator()
	{
		if (VERBOSE > 0) printf("Destructor of operator term!");
	}

	void applyProx(flexBoxData<T, Tvector>* data, const std::vector<T> &sigma, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		switch (this->p)
		{
			case dualL1AnisoProx:
			{
				flexProxList<T, Tvector>::dualL1AnisoProx(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
			case dualL1IsoProx:
			{
				flexProxList<T, Tvector>::dualL1IsoProx(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
			case dualL2Prox:
			{
				flexProxList<T, Tvector>::dualL2Prox(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
			case dualFrobeniusProx:
			{
				flexProxList<T, Tvector>::dualFrobeniusProx(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
		}
	};

#if __CUDACC__
	__device__ void applyProxElement(T** ptrYList, const T* yTilde, T** sigmaList, const int* dualNumbers, int numDuals, int index)
	{
		switch (this->p)
		{
			case dualL1AnisoProx:
			{
				T alphaTmp = this->alpha;

				for (int i = 0; i < numDuals; ++i)
				{
					const int dualNum = dualNumbers[i];

					ptrYList[dualNum][index] = myMinGPUf(alphaTmp, myMaxGPUf(-alphaTmp, yTilde[dualNum]));
				}
				break;
			}
			case dualL1IsoProx:
			{
				//square yTilde
				T sumTildeSquared = (T)0;

				for (int i = 0; i < numDuals; ++i)
				{
					const int dualNum = dualNumbers[i]; 
					
					sumTildeSquared += yTilde[dualNum] * yTilde[dualNum];
				}

				sumTildeSquared = (T)1 / myMaxGPUf((T)1, sqrtf(sumTildeSquared) / this->alpha);

				for (int i = 0; i < numDuals; ++i)
				{
					const int dualNum = dualNumbers[i];

					ptrYList[dualNum][index] = yTilde[dualNum] * sumTildeSquared;
				}

				break;
			}
			case dualHuberProx:
			{
				//square yTilde
				T sumTildeSquared = (T)0;
				T epsi = (T)0.01;

				for (int i = 0; i < numDuals; ++i)
				{
					const int dualNum = dualNumbers[i];

					T tmpyTilde = yTilde[dualNum] / ((T)1 + sigmaList[dualNum][index] * (T)0.01 / this->alpha);
					sumTildeSquared += tmpyTilde * tmpyTilde;
				}

				sumTildeSquared = (T)1 / myMaxGPUf((T)1, sqrtf(sumTildeSquared) / this->alpha);

				for (int i = 0; i < numDuals; ++i)
				{
					const int dualNum = dualNumbers[i];

					ptrYList[dualNum][index] = yTilde[dualNum] * sumTildeSquared;
				}

				break;
			}
			case dualL2Prox:
			{
				for (int i = 0; i < numDuals; ++i)
				{
					const int dualNum = dualNumbers[i];

					ptrYList[dualNum][index] = this->alpha / (sigmaList[dualNum][index] + this->alpha) * yTilde[dualNum];
				}
				break;
			}

		}
	}
#endif
};

#endif