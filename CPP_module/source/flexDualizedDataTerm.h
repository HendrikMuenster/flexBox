#ifndef flexDualizedDataTerm_H
#define flexDualizedDataTerm_H

#include "flexTermDual.h"
#include "prox/flexProx.h"

template < typename T, typename Tvector >
class flexDualizedDataTerm : public flexTermDual<T, Tvector>
{
private:

public:
	Tvector f;
	std::vector<std::vector<T>> fList;
	flexProx<T,Tvector>* myProx;
	
	#if __CUDACC__
		//save pointer to device memory
		T* fPtr;
	#endif
	
		flexDualizedDataTerm(flexProx<T,Tvector>* _myProx, T _alpha, int numberPrimals, std::vector<flexLinearOperator<T, Tvector>* > _operatorList) : flexDualizedDataTerm(_myProx,_alpha,numberPrimals,_operatorList, std::vector<std::vector<T>>(0)){};

		flexDualizedDataTerm(flexProx<T,Tvector>* _myProx, T _alpha, int numberPrimals, std::vector<flexLinearOperator<T, Tvector>* > _operatorList, std::vector<std::vector<T>> _fList) : myProx(_myProx), flexTermDual<T, Tvector>(_myProx->getProx(), _alpha, (int)_operatorList.size(), (int)_operatorList.size() / numberPrimals)	
		{
			//convert input f to device vector
			fList.resize(_fList.size());

			#if __CUDACC__
                //todo! f is now list of data. Implement for CUDA
            
                this->f.resize(_fList[0].size());
				thrust::copy(_fList[0].begin(), _fList[0].end(), this->f.begin());
				fPtr = thrust::raw_pointer_cast(this->f.data());
			#else
				for (int i = 0; i < fList.size(); ++i)
				{
					this->fList[i].resize(_fList[i].size());
					
					std::copy(_fList[i].begin(), _fList[i].end(), this->fList[i].begin());
				}
			#endif

			this->operatorList = _operatorList;

			//create sigma and tau
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

	~flexDualizedDataTerm()
	{
		delete myProx;
		if (VERBOSE > 0) printf("Destructor of data term!");
	}

	void applyProx(flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
#if __CUDACC__
#else
		
		switch (this->p)
		{
			case dualL2DataProx:
			case dualL1DataProx:
			case dualKLDataProx:
			{
				myProx->applyProx(this->alpha, data,dualNumbers,primalNumbers,this->fList);
				break;
			}
			case dualL1AnisoProx:
			case dualL1IsoProx:
			case dualL2Prox:
			case dualFrobeniusProx:
			case dualBoxConstraintProx:
			{
				myProx->applyProx(this->alpha, data,dualNumbers,primalNumbers);
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
			case dualL2DataProx:
			{
				for (int i = 0; i < numDuals; ++i)
				{
					const int dualNum = dualNumbers[i];

					T tmpSigma = sigmaList[dualNum][index];

					T factor = this->alpha / (tmpSigma + this->alpha);
					ptrYList[dualNum][index] = factor * (yTilde[dualNum] - tmpSigma * this->fPtr[index]);
				}
				break;
			}
			case dualL1DataProx:
			{
				for (int i = 0; i < numDuals; ++i)
				{
					const int dualNum = dualNumbers[i];

					ptrYList[dualNum][index] = myMinGPU(this->alpha, myMaxGPU(-this->alpha, yTilde[dualNum] - this->alpha * sigmaList[dualNum][index] * this->fPtr[index]));
				}
				break;
			}
			case dualKLDataProx:
			{
				
				break;
			}
		}
	}
#endif
};

#endif