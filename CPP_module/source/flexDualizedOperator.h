#ifndef flexDualizedOperator_H
#define flexDualizedOperator_H

#include "flexTermDual.h"

template < typename T, typename Tvector >
class flexDualizedOperator : public flexTermDual<T, Tvector>
{
public:
	flexDualizedOperator(prox _p, T _alpha, int numberPrimals, std::vector<flexLinearOperator<T, Tvector>* > _operatorList) : flexTermDual<T, Tvector>(_p, _alpha, numberPrimals, _operatorList.size() / numberPrimals, 0.0f)
	{
		this->operatorList = _operatorList;

		for (int i = 0; i < _operatorList.size() / numberPrimals; ++i)
		{
			for (int j = 0; j < numberPrimals; ++j)
			{
				int opNum = i*numberPrimals + j;

				this->operatorListT.push_back(this->operatorList[opNum]->copy());
				this->operatorListT[opNum]->transpose();
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

	void applyProx(flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
#if __CUDACC__
#else
		switch (this->p)
		{
			case dualL1AnisoProx:
			{
				for (int k = 0; k < dualNumbers.size(); k++)
				{
					T* ptrY = data->y[dualNumbers[k]].data();
					T* ptrYtilde = data->yTilde[dualNumbers[k]].data();

					int numElements = data->yTilde[dualNumbers[k]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						ptrY[i] = myMin<T>(this->alpha, myMax<T>(-this->alpha, ptrYtilde[i]));
					}
				}

				//flexProxList<T, Tvector>::dualL1AnisoProx(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
			case dualL1IsoProx:
			{
				if (dualNumbers.size() == 1)
				{
					T* ptrY0 = data->y[dualNumbers[0]].data();

					T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();

					int numElements = data->yTilde[dualNumbers[0]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						T yTmp = (T)1 / myMax<T>((T)1, fabs(ptrYtilde0[i]) / this->alpha);

						ptrY0[i] = ptrYtilde0[i] * yTmp;
					}
				}
				else if (dualNumbers.size() == 2)
				{
					T* ptrY0 = data->y[dualNumbers[0]].data();
					T* ptrY1 = data->y[dualNumbers[1]].data();

					T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
					T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();

					int numElements = data->yTilde[dualNumbers[0]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i]) + pow2(ptrYtilde1[i])) / this->alpha);

						ptrY0[i] = ptrYtilde0[i] * yTmp;
						ptrY1[i] = ptrYtilde1[i] * yTmp;
					}
				}
				else if (dualNumbers.size() == 3)
				{
					T* ptrY0 = data->y[dualNumbers[0]].data();
					T* ptrY1 = data->y[dualNumbers[1]].data();
					T* ptrY2 = data->y[dualNumbers[2]].data();

					T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
					T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();
					T* ptrYtilde2 = data->yTilde[dualNumbers[2]].data();

					int numElements = data->yTilde[dualNumbers[0]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i]) + pow2(ptrYtilde1[i]) + pow2(ptrYtilde2[i])) / this->alpha);

						ptrY0[i] = ptrYtilde0[i] * yTmp;
						ptrY1[i] = ptrYtilde1[i] * yTmp;
						ptrY2[i] = ptrYtilde2[i] * yTmp;
					}
				}
				else
				{
					printf("Alert! Iso prox not implemented for dim>3");
				}

				//flexProxList<T, Tvector>::dualL1IsoProx(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
			case dualHuberProx:
			{
				if (dualNumbers.size() == 1)
				{
					T* ptrY0 = data->y[dualNumbers[0]].data();

					T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();

					T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

					int numElements = data->yTilde[dualNumbers[0]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * (T)0.01 / this->alpha);

						T yTmp = (T)1 / myMax<T>((T)1, fabs(ptrYtilde0[i] * huberFactor) / this->alpha);

						ptrY0[i] = ptrYtilde0[i] * yTmp;
					}
				}
				else if (dualNumbers.size() == 2)
				{
					T* ptrY0 = data->y[dualNumbers[0]].data();
					T* ptrY1 = data->y[dualNumbers[1]].data();

					T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
					T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();

					T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

					int numElements = data->yTilde[dualNumbers[0]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * (T)0.01 / this->alpha);

						T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i] * huberFactor) + pow2(ptrYtilde1[i] * huberFactor)) / this->alpha);

						ptrY0[i] = ptrYtilde0[i] * yTmp;
						ptrY1[i] = ptrYtilde1[i] * yTmp;
					}
				}
				else if (dualNumbers.size() == 3)
				{
					T* ptrY0 = data->y[dualNumbers[0]].data();
					T* ptrY1 = data->y[dualNumbers[1]].data();
					T* ptrY2 = data->y[dualNumbers[2]].data();

					T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
					T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();
					T* ptrYtilde2 = data->yTilde[dualNumbers[2]].data();

					T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

					int numElements = data->yTilde[dualNumbers[0]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * (T)0.01 / this->alpha);

						T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i] * huberFactor) + pow2(ptrYtilde1[i] * huberFactor) + pow2(ptrYtilde2[i] * huberFactor)) / this->alpha);

						ptrY0[i] = ptrYtilde0[i] * yTmp;
						ptrY1[i] = ptrYtilde1[i] * yTmp;
						ptrY2[i] = ptrYtilde2[i] * yTmp;
					}
				}
				else
				{
					printf("Alert! Huber prox not implemented for dim>3");
				}

				//flexProxList<T, Tvector>::dualHuberProx(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
			case dualL2Prox:
			{
				for (int k = 0; k < dualNumbers.size(); k++)
				{
					T* ptrY = data->y[dualNumbers[k]].data();
					T* ptrYtilde = data->yTilde[dualNumbers[k]].data();

					T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

					int numElements = data->yTilde[dualNumbers[k]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						ptrY[i] = this->alpha / (ptrSigma[i] + this->alpha) * ptrYtilde[i];
					}
				}

				//flexProxList<T, Tvector>::dualL2Prox(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
			case dualFrobeniusProx:
			{
				T norm = (T)0;
				for (int k = 0; k < dualNumbers.size(); k++)
				{
					T* ptrYTilde = data->yTilde[dualNumbers[k]].data();

					int numElements = data->yTilde[dualNumbers[k]].size();

					#pragma omp parallel for reduction(+: norm)
					for (int i = 0; i < numElements; i++)
					{
						norm += ptrYTilde[i] * ptrYTilde[i];
					}
				}

				norm = (T)1 / std::max((T)1, std::sqrt(norm) / this->alpha);

				for (int k = 0; k < dualNumbers.size(); k++)
				{
					T* ptrY = data->y[dualNumbers[k]].data();
					T* ptrYTilde = data->yTilde[dualNumbers[k]].data();

					int numElements = data->yTilde[dualNumbers[k]].size();

					#pragma omp parallel for
					for (int i = 0; i < numElements; i++)
					{
						ptrY[i] = ptrYTilde[i] * norm;
					}
				}

				//flexProxList<T, Tvector>::dualFrobeniusProx(data, sigma, dualNumbers, primalNumbers, this->alpha);
				break;
			}
		}
#endif
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