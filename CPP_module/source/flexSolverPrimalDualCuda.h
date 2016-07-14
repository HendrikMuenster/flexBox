#ifndef flexSolverPrimalDualCuda_H
#define flexSolverPrimalDualCuda_H

#include "flexSolver.h"

#include "kernelPrimalDual.h"


template < class T, class Tvector >
class flexSolverPrimalDualCuda : public  flexSolver<T,Tvector>
{
private:
	//Params
	T theta;

	//List of primal terms is list of pointers to terms
	std::vector<flexTermPrimal<T, Tvector>*> termsPrimal;
	//List of dual terms is list of pointers to terms
	std::vector<flexTermDual<T, Tvector>*> termsDual;

	//create GPU version of everything
	thrust::device_vector<flexTermPrimal<T, Tvector>*> termsPrimalG;
	thrust::device_vector<flexTermDual<T, Tvector>*> termsDualG;
	thrust::device_vector<int* > pcpG;
	thrust::device_vector<int* > dcpG;
	thrust::device_vector<int* > dcdG;

	thrust::device_vector<T*> listX, listXOld, listXTilde, listXError, listXBar;
	thrust::device_vector<T*> listY, listYOld, listYTilde, listYError;

	thrust::device_vector<T*> listTauElt;
	thrust::device_vector<T*> listSigmaElt;

	thrust::device_vector<int> listSizeX, listSizeY;

	T** listXptr;
	T** listXOldptr;
	T** listXErrorptr;
	T** listXBarptr;
	T** listXTildeptr;

	T** listYptr;
	T** listYOldptr;
	T** listYErrorptr;
	T** listYTildeptr;

	int gridSizeX;
	int gridSizeY;
	int blockSize;

public:

	flexSolverPrimalDualCuda()
	{
		this->theta = static_cast<T>(1);
	}

	~flexSolverPrimalDualCuda()
	{
		if (VERBOSE > 0) printf("Destructor flexSolverPrimalDualCuda solver\n!");

		for (int i = termsPrimal.size() - 1; i >= 0; --i)
		{
			cudaFree(termsPrimalG[i]); //CUDA_CHECK;
			delete termsPrimal[i];
		}

		for (int i = termsDual.size() - 1; i >= 0; --i)
		{
			cudaFree(termsDualG[i]); //CUDA_CHECK;
			delete termsDual[i];
		}

		for (int i = pcpG.size() - 1; i >= 0; --i)
		{
			cudaFree(pcpG[i]); //CUDA_CHECK;
		}

		for (int i = dcpG.size() - 1; i >= 0; --i)
		{
			cudaFree(dcpG[i]); //CUDA_CHECK;
			cudaFree(dcdG[i]); //CUDA_CHECK;
		}

		listX.clear();
		listXOld.clear();
		listXError.clear();
		listXBar.clear();
		listY.clear();
		listYOld.clear();
		listYError.clear();
	}

	void init(flexBoxData<T, Tvector> *data)
	{
		//zero init
		int maxSizeX = 0;
		for (int i = 0; i < data->x.size(); ++i)
		{
			listX.push_back(gpuPtr(data->x[i]));
			listXOld.push_back(gpuPtr(data->xOld[i]));
			listXError.push_back(gpuPtr(data->xError[i]));
			listXBar.push_back(gpuPtr(data->xBar[i]));
			listXTilde.push_back(gpuPtr(data->xTilde[i]));

			listTauElt.push_back(gpuPtr(data->tauElt[i]));
			listSizeX.push_back(data->x[i].size());

			maxSizeX = myMax<int>(maxSizeX, data->x[i].size());
		}


		listXptr = gpuPtr(listX);
		listXOldptr = gpuPtr(listXOld);
		listXErrorptr = gpuPtr(listXError);
		listXBarptr = gpuPtr(listXBar);
		listXTildeptr = gpuPtr(listXTilde);

		int maxSizeY = 0;
		for (int i = 0; i < data->y.size(); ++i)
		{
			listY.push_back(gpuPtr(data->y[i]));
			listYOld.push_back(gpuPtr(data->yOld[i]));
			listYError.push_back(gpuPtr(data->yError[i]));
			listYTilde.push_back(gpuPtr(data->yTilde[i]));

			listSigmaElt.push_back(gpuPtr(data->sigmaElt[i]));
			listSizeY.push_back(data->y[i].size());

			maxSizeY = myMax<int>(maxSizeY, data->y[i].size());

		}

		listYptr = gpuPtr(listY);
		listYOldptr = gpuPtr(listYOld);
		listYErrorptr = gpuPtr(listYError);
		listYTildeptr = gpuPtr(listYTilde);

		this->blockSize = 512;

		this->gridSizeX = (maxSizeX + this->blockSize - 1) / this->blockSize;
		this->gridSizeY = (maxSizeY + this->blockSize - 1) / this->blockSize;
		//calculate element-wise tau and sigma
		
		tauElement << <this->gridSizeX, this->blockSize >> >(termsDual.size(), data->x.size(), gpuPtr(termsDualG), gpuPtr(listTauElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(listSizeX)); //CUDA_CHECK;
		sigmaElement << <this->gridSizeY, this->blockSize >> >(termsDual.size(), data->y.size(), gpuPtr(termsDualG), gpuPtr(listSigmaElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(listSizeY)); //CUDA_CHECK;



	}

	void addPrimal(flexTermPrimal<T, Tvector>* _primalPart, std::vector<int> _correspondingPrimals)
	{
		termsPrimal.push_back(_primalPart); 
		
		flexTermPrimal<T, Tvector>* ptrPrimal;
		switch (_primalPart->p)
		{
			case primalEmptyProx:
			{
				cudaMalloc((void **)&ptrPrimal, sizeof(flexTermPrimal<T, Tvector>)); //CUDA_CHECK;
				cudaMemcpy(ptrPrimal, _primalPart, sizeof(flexTermPrimal<T, Tvector>), cudaMemcpyHostToDevice);
				break;
			}
		}
		termsPrimalG.push_back(ptrPrimal);
		pcpG.push_back(gpuPtr(new thrust::device_vector<int>(_correspondingPrimals)));
	}

	void addDual(flexBoxData<T, Tvector> *data,flexTermDual<T, Tvector>* _dualPart, std::vector<int> _correspondingPrimals)
	{
		termsDual.push_back(_dualPart);

		std::vector<int> tmpDCD;
		for (int i = 0; i < _dualPart->getNumberVars(); ++i)
		{
			data->addDualVar(_dualPart->dualVarLength(i));

			tmpDCD.push_back(data->getNumDualVars() - 1);
		}
			
		flexTermDual<T, Tvector>* ptrDual;
		switch (_dualPart->p)
		{
		case dualL1AnisoProx:
		case dualL1IsoProx:
		case dualL2Prox:
		case dualFrobeniusProx:
		case dualHuberProx:
		{
			cudaMalloc((void **)&ptrDual, sizeof(flexDualizedOperator<T, Tvector>)); //CUDA_CHECK;
			cudaMemcpy(ptrDual, _dualPart, sizeof(flexDualizedOperator<T, Tvector>), cudaMemcpyHostToDevice);
			break;
		}
		case dualL2DataProx:
		case dualL1DataProx:
		case dualKLDataProx:
		{
			cudaMalloc((void **)&ptrDual, sizeof(flexDualizedDataTerm<T, Tvector>)); //CUDA_CHECK;
			cudaMemcpy(ptrDual, _dualPart, sizeof(flexDualizedDataTerm<T, Tvector>), cudaMemcpyHostToDevice);
			break;
		}
		}
		termsDualG.push_back(ptrDual);


		dcpG.push_back(gpuPtr(new thrust::device_vector<int>(_correspondingPrimals)));
		dcdG.push_back(gpuPtr(new thrust::device_vector<int>(tmpDCD)));
	}

	void doIteration(flexBoxData<T, Tvector> *data)
	{
		std::swap(listXptr, listXOldptr);
		std::swap(listYptr, listYOldptr);

		updateY << <this->gridSizeY, this->blockSize >> >(termsDual.size(), data->y.size(), gpuPtr(termsDualG), listYptr, listYOldptr, listXBarptr, gpuPtr(listSigmaElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(listSizeY)); //CUDA_CHECK;
		
		/*switch (data->x.size())
		{
			case 1:
			{
				updateX<T, Tvector,1 > << <this->gridSizeX, this->blockSize >> >(termsDual.size(), termsPrimal.size(), data->x.size(), gpuPtr(termsDualG), gpuPtr(termsPrimalG), listXptr, listXOldptr, listXTildeptr, listXBarptr, listYptr, gpuPtr(listTauElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(pcpG), gpuPtr(listSizeX)); CUDA_CHECK;
				break;
			}
			case 2:
			{
				updateX<T, Tvector,2 > << <this->gridSizeX, this->blockSize >> >(termsDual.size(), termsPrimal.size(), data->x.size(), gpuPtr(termsDualG), gpuPtr(termsPrimalG), listXptr, listXOldptr, listXTildeptr, listXBarptr, listYptr, gpuPtr(listTauElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(pcpG), gpuPtr(listSizeX)); CUDA_CHECK;
				break;
			}
			case 3:
			{
				updateX<T, Tvector,3 > << <this->gridSizeX, this->blockSize >> >(termsDual.size(), termsPrimal.size(), data->x.size(), gpuPtr(termsDualG), gpuPtr(termsPrimalG), listXptr, listXOldptr, listXTildeptr, listXBarptr, listYptr, gpuPtr(listTauElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(pcpG), gpuPtr(listSizeX)); CUDA_CHECK;
				break;
			}
		}*/

		updateX<< <this->gridSizeX, this->blockSize >> >(termsDual.size(), termsPrimal.size(), data->x.size(), gpuPtr(termsDualG), gpuPtr(termsPrimalG), listXptr, listXOldptr, listXTildeptr, listXBarptr, listYptr, gpuPtr(listTauElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(pcpG), gpuPtr(listSizeX)); //CUDA_CHECK;

	}

	T calculateError(flexBoxData<T, Tvector> *data)
	{
		errorY << <this->gridSizeY, this->blockSize >> >(termsDual.size(), data->y.size(), gpuPtr(termsDualG), listYptr, listYOldptr, listYErrorptr, listXptr, listXOldptr, gpuPtr(listSigmaElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(listSizeY)); //CUDA_CHECK;
		errorX << <this->gridSizeX, this->blockSize >> >(termsDual.size(), data->y.size(), data->x.size(), gpuPtr(termsDualG), listYptr, listYOldptr, listXErrorptr, listXptr, listXOldptr, gpuPtr(listTauElt), gpuPtr(dcdG), gpuPtr(dcpG), gpuPtr(listSizeX)); //CUDA_CHECK;

		//sum things up
		T primalResidual = static_cast<T>(0);
		T dualResidual = static_cast<T>(0);

		for (int i = 0; i < data->getNumPrimalVars(); ++i)
		{
			vectorAbs(data->xError[i]);

			primalResidual += vectorSum(data->xError[i]) / static_cast<T>(data->xError[i].size());

		}
		primalResidual = primalResidual / (T)data->getNumPrimalVars();

		for (int i = 0; i < data->getNumDualVars(); ++i)
		{
			vectorAbs(data->yError[i]); //take absolute value of entries

			dualResidual += vectorSum(data->yError[i]) / static_cast<T>(data->yError[i].size()); //sum values up and add to dual residual
		}
		dualResidual = dualResidual / (T)data->getNumDualVars();

		return primalResidual + dualResidual;
	}

	T* gpuPtr(Tvector &input)
	{
		return thrust::raw_pointer_cast(input.data());
	}

	int* gpuPtr(thrust::device_vector<int> &input)
	{
		return thrust::raw_pointer_cast(input.data());
	}

	int* gpuPtr(thrust::device_vector<int>* input)
	{
		return thrust::raw_pointer_cast(input->data());
	}

	int** gpuPtr(thrust::device_vector<int*> &input)
	{
		return thrust::raw_pointer_cast(input.data());
	}

	T** gpuPtr(thrust::device_vector<T*> &input)
	{
		return thrust::raw_pointer_cast(input.data());
	}

	flexLinearOperator < T, Tvector >** gpuPtr(thrust::device_vector<flexLinearOperator < T, Tvector >*> &input)
	{
		return thrust::raw_pointer_cast(input.data());
	}

	linOp* gpuPtr(thrust::device_vector<linOp> &input)
	{
		return thrust::raw_pointer_cast(input.data());
	}

	flexTermDual<T, Tvector>** gpuPtr(thrust::device_vector<flexTermDual<T, Tvector>*> &input)
	{
		return thrust::raw_pointer_cast(input.data());
	}

	flexTermPrimal<T, Tvector>** gpuPtr(thrust::device_vector<flexTermPrimal<T, Tvector>*> &input)
	{
		return thrust::raw_pointer_cast(input.data());
	}
};

#endif