#ifndef flexSolverPrimalDual_H
#define flexSolverPrimalDual_H

#include "flexSolver.h"

template < class T, class Tvector >
class flexSolverPrimalDual : public  flexSolver<T, Tvector>
{
private:
	//Params
	T theta;

	//List of primal terms is list of pointers to terms
	std::vector<flexTermPrimal<T, Tvector>*> termsPrimal;
	//List of dual terms is list of pointers to terms
	std::vector<flexTermDual<T, Tvector>*> termsDual;

	//List of primal variables corresponding to primal terms
	std::vector<std::vector<int> > pcp;
	//List of primal variables corresponding to dual terms
	std::vector<std::vector<int> > dcp;
	//List of dual variables corresponding to dual terms
	std::vector<std::vector<int> > dcd;
public:

	flexSolverPrimalDual()
	{
		this->theta = static_cast<T>(1);
	}

	~flexSolverPrimalDual()
	{
		for (int i = termsPrimal.size() - 1; i >= 0; --i)
		{
			delete termsPrimal[i];
		}

		for (int i = termsDual.size() - 1; i >= 0; --i)
		{
			delete termsDual[i];
		}
	}

	void init(flexBoxData<T, Tvector> *data)
	{
		this->calculateTauSigma(data);
	}

	void calculateTauSigma(flexBoxData<T, Tvector> *data)
	{
		std::vector<T> tmpVec;

		for (int k = 0; k < termsDual.size(); ++k)
		{
			auto numPrimals = dcp[k].size();
			auto numDuals = dcd[k].size();

			for (int i = 0; i < numDuals; ++i)
			{
				const int dualNum = dcd[k][i];

				for (int j = 0; j < numPrimals; ++j)
				{
					const int primalNum = dcp[k][j];

					auto operatorNumber = numPrimals * i + j;

					tmpVec = termsDual[k]->operatorList[operatorNumber]->getAbsRowSum();
					vectorPlus(data->sigmaElt[dualNum], tmpVec);

					tmpVec = termsDual[k]->operatorListT[operatorNumber]->getAbsRowSum();
					vectorPlus(data->tauElt[primalNum], tmpVec);
				}
			}
		}

		for (int i = 0; i < data->y.size(); ++i)
		{
			auto ptrSigma = data->sigmaElt[i].data();

			int numElements = data->y[i].size();

			#pragma omp parallel for
			for (int j = 0; j < numElements; ++j)
			{
				ptrSigma[j] = (T)1 / ptrSigma[j];
			}
		}

		for (int i = 0; i < data->x.size(); ++i)
		{
			auto ptrTau = data->tauElt[i].data();

			int numElements = data->x[i].size();

			#pragma omp parallel for
			for (int j = 0; j < numElements; ++j)
			{
				ptrTau[j] = (T)1 / ptrTau[j];
			}
		}
	}

	void addPrimal(flexTermPrimal<T, Tvector>* _primalPart, std::vector<int> _correspondingPrimals)
	{
		this->termsPrimal.push_back(_primalPart);

		this->pcp.push_back(_correspondingPrimals);
	}

	void addDual(flexBoxData<T, Tvector> *data, flexTermDual<T, Tvector>* _dualPart, std::vector<int> _correspondingPrimals)
	{
		termsDual.push_back(_dualPart);

		std::vector<int> tmpDCD;
		for (int i = 0; i < _dualPart->getNumberVars(); ++i)
		{
			data->addDualVar(_dualPart->dualVarLength(i));

			tmpDCD.push_back(data->getNumDualVars() - 1);
		}

		dcp.push_back(_correspondingPrimals);
		dcd.push_back(tmpDCD);
	}

	void yTilde(flexBoxData<T, Tvector>* data, flexTermDual<T, Tvector>* dualTerm, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
            const int dualNum = dualNumbers[i];
            
            vectorScalarSet(data->yTilde[dualNum],(T)0);
            
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				int operatorNumber = primalNumbers.size() * i + j;

				dualTerm->operatorList[operatorNumber]->timesPlus(data->xBar[primalNumbers[j]], data->yTilde[dualNum]);
			}
            
            T* ptrYtilde = data->yTilde[dualNum].data();
            T* ptrYold = data->yOld[dualNum].data();
			T* ptrSigma = data->sigmaElt[dualNum].data();
            
            int numElements = data->yTilde[dualNum].size();
            
            #pragma omp parallel for
            for (int k = 0; k < numElements; ++k)
            {
				ptrYtilde[k] = ptrYold[k] + ptrSigma[k] * ptrYtilde[k];
            }
		}
	}

	void xTilde(flexBoxData<T, Tvector>* data, flexTermDual<T, Tvector>* dualTerm, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				int operatorNumber = primalNumbers.size() * i + j;

				const int primalNum = primalNumbers[j];
				const int dualNum = dualNumbers[i];

				dualTerm->operatorListT[operatorNumber]->timesPlus(data->y[dualNum], data->xTilde[primalNum]);
			}
		}
	}

	void doIteration(flexBoxData<T, Tvector> *data)
	{
		//bool doTime = false;

		//Timer timer;

		//timer.reset();
		data->xOld.swap(data->x);
		data->yOld.swap(data->y);

		//timer.end(); printf("Time for swap was: %f\n", timer.elapsed());

		//if (doTime) timer.reset();
		for (int i = 0; i < termsDual.size(); ++i)
		{
			//timer.reset();
			this->yTilde(data, termsDual[i], dcd[i], dcp[i]);
			//timer.end(); printf("Time for termsDual[i]->yTilde(data,sigma, dcd[i], dcp[i]); was: %f\n", timer.elapsed());
			
			//timer.reset();
			termsDual[i]->applyProx(data, dcd[i], dcp[i]);
			//timer.end();printf("Time for termsDual[i]->applyProx(data,sigma, dcd[i], dcp[i]); was: %f\n", timer.elapsed());
		}

		//timer.reset();
		for (int i = 0; i < data->xTilde.size(); ++i)
		{
			vectorScalarSet(data->xTilde[i], (T)0);
		}
		//timer.end(); printf("Time for vectorScalarSet(data->xTilde[i], (T)0); was: %f\n", timer.elapsed());

		//timer.reset();
		for (int i = 0; i < termsDual.size(); ++i)
		{
			this->xTilde(data, termsDual[i], dcd[i], dcp[i]);
		}
		// set tildeX = xOld - tau * tildeX
		for (int i = 0; i < data->xTilde.size(); ++i)
		{
			T* ptrXtilde = data->xTilde[i].data();
			T* ptrXold = data->xOld[i].data();
			T* ptrTau = data->tauElt[i].data();

			int numElements = data->xTilde[i].size();

			#pragma omp parallel for
			for (int k = 0; k < numElements; ++k)
			{
				ptrXtilde[k] = ptrXold[k] - ptrTau[k] * ptrXtilde[k];
			}
		}
		//timer.end(); printf("Time for this->xTilde was: %f\n", timer.elapsed());

		//timer.reset();
		for (int i = 0; i < termsPrimal.size(); ++i)
		{
			termsPrimal[i]->applyProx(data, pcp[i]);
		}
		//timer.end(); printf("Time for termsPrimal[i]->applyProx(data, tau, pcp[i]); was: %f\n", timer.elapsed());

		//do overrelaxation
		//timer.reset();
		for (int i = 0; i < data->xTilde.size(); ++i)
		{
			T* ptrXbar = data->xBar[i].data();
			T* ptrX = data->x[i].data();
			T* ptrXold = data->xOld[i].data(); 
			
			int numElements = data->x[i].size();

			#pragma omp parallel for
			for (int k = 0; k < numElements; ++k)
			{
				ptrXbar[k] = 2 * ptrX[k] - ptrXold[k];
			}
		}
		//timer.end(); printf("Time for doOverrelaxation(data->x[i], data->xOld[i], data->xBar[i]); was: %f\n", timer.elapsed());
	}

	void yError(flexBoxData<T, Tvector>* data, flexTermDual<T, Tvector>* dualTerm, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				//Check here
				int operatorNumber = primalNumbers.size() * i + j;

				const int primalNum = primalNumbers[j];
				const int dualNum = dualNumbers[i];

				data->xTmp[primalNum] = data->x[primalNum];
				vectorMinus(data->xTmp[primalNum], data->xOld[primalNum]);

				dualTerm->operatorList[operatorNumber]->timesMinus(data->xTmp[primalNum], data->yError[dualNum]);
			}
		}
	}

	void xError(flexBoxData<T, Tvector>* data, flexTermDual<T, Tvector>* dualTerm, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		for (int i = 0; i < dualNumbers.size(); ++i)
		{
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				//Check here
				int operatorNumber = primalNumbers.size() * i + j;

				const int primalNum = primalNumbers[j];
				const int dualNum = dualNumbers[i];

				data->yTmp[dualNum] = data->y[dualNum];
				vectorMinus(data->yTmp[dualNum], data->yOld[dualNum]);

				dualTerm->operatorListT[operatorNumber]->timesMinus(data->yTmp[dualNum], data->xError[primalNum]);
			}
		}
	}

	T calculateError(flexBoxData<T, Tvector> *data)
	{
		//first part of primal and dual residuals
		for (int i = 0; i < data->x.size(); ++i)
		{
			T* ptrxError = data->xError[i].data();
			T* ptrX = data->x[i].data();
			T* ptrXold = data->xOld[i].data();

			T* ptrTau = data->tauElt[i].data();

			int numElements = data->x[i].size();

			#pragma omp parallel for
			for (int k = 0; k < numElements; ++k)
			{
				ptrxError[k] = (ptrX[k] - ptrXold[k]) / ptrTau[k];
			}
		}

		for (int i = 0; i < data->y.size(); ++i)
		{
			T* ptrYError = data->yError[i].data();
			T* ptrY = data->y[i].data();
			T* ptrYold = data->yOld[i].data();

			T* ptrSigma = data->sigmaElt[i].data();

			int numElements = data->y[i].size();

			#pragma omp parallel for
			for (int k = 0; k < numElements; ++k)
			{
				ptrYError[k] = (ptrY[k] - ptrYold[k]) / ptrSigma[k];
			}
		}

		//operator specifiy error
		for (int i = 0; i < termsDual.size(); ++i)
		{
			this->xError(data, termsDual[i], dcd[i], dcp[i]);
			this->yError(data, termsDual[i], dcd[i], dcp[i]);
		}

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

};

#endif