#ifndef flexSolverPrimalDual_H
#define flexSolverPrimalDual_H

#include "flexSolver.h"

template < class T, class Tvector >
class flexSolverPrimalDual : public  flexSolver<T, Tvector>
{
private:
	//Params
	T theta;

	std::vector<T> tau;
	std::vector<T> sigma;

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
		//zero init
		for (int i = 0; i < data->x.size(); ++i)
		{
			tau.push_back(static_cast<T>(0));
		}

		for (int i = 0; i < data->y.size(); ++i)
		{
			sigma.push_back(static_cast<T>(0));
		}

		for (int i = 0; i < termsDual.size(); ++i)
		{
			for (int j = 0; j < dcp[i].size(); ++j)
			{
				tau[dcp[i][j]] += termsDual[i]->myTau[j];
			}

			for (int j = 0; j < dcd[i].size(); ++j)
			{
				sigma[dcd[i][j]] += termsDual[i]->mySigma[j];
			}
		}

		for (int i = 0; i < tau.size(); ++i)
		{
			tau[i] = static_cast<T>(1) / tau[i];
		}

		for (int i = 0; i < sigma.size(); ++i)
		{
			sigma[i] = static_cast<T>(1) / sigma[i];
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
			for (int j = 0; j < primalNumbers.size(); ++j)
			{
				int operatorNumber = primalNumbers.size() * i + j;

				//Tvector xBarTmp(data->xBar[primalNumbers[j]]);
				//vectorScalarProduct(xBarTmp, sigma[dualNumbers[i]]);

				const int primalNum = primalNumbers[j];
				const int dualNum = dualNumbers[i];

				data->xTmp[primalNum] = data->xBar[primalNum];
				vectorScalarProduct(data->xTmp[primalNum], sigma[dualNum]);

				dualTerm->operatorList[operatorNumber]->timesPlus(data->xTmp[primalNum], data->yTilde[dualNum]);
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

				data->yTmp[dualNum] = data->y[dualNum];
				vectorScalarProduct(data->yTmp[dualNum], tau[primalNum]);

				dualTerm->operatorListT[operatorNumber]->timesMinus(data->yTmp[dualNum], data->xTilde[primalNum]);
			}
		}
	}

	void doIteration(flexBoxData<T, Tvector> *data)
	{
		//bool doTime = false;

		//Timer timer;

		//if (doTime) timer.reset();
		data->xOld.swap(data->x);
		data->yOld.swap(data->y);

		//if (doTime) { timer.end(); printf("Time for swap was: %f\n", timer.elapsed()); }

		//if (doTime) timer.reset();
		for (int i = 0; i < data->yTilde.size(); ++i)
		{
			data->yTilde[i] = data->yOld[i];
		}


		for (int i = 0; i < termsDual.size(); ++i)
		{
			this->yTilde(data, termsDual[i], dcd[i], dcp[i]);
			//if (doTime) timer.reset();
			//termsDual[i]->yTilde(data, sigma, dcd[i], dcp[i]);
			//if (doTime) timer.end(); 
			//if (doTime) printf("Time for termsDual[i]->yTilde(data,sigma, dcd[i], dcp[i]); was: %f\n", timer.elapsed());
			//if (doTime) timer.reset();
			termsDual[i]->applyProx(data, sigma, dcd[i], dcp[i]);
			//if (doTime) timer.end();
			//if (doTime) printf("Time for termsDual[i]->applyProx(data,sigma, dcd[i], dcp[i]); was: %f\n", timer.elapsed());

		}
		//if (doTime) timer.reset();
		for (int i = 0; i < data->xTilde.size(); ++i)
		{
			data->xTilde[i] = data->xOld[i];
		}
		//if (doTime) timer.end();
		//if (doTime) printf("Time for data->xTilde[i] = data->xOld[i]; was: %f\n", timer.elapsed());

		//if (doTime) timer.reset();
		for (int i = 0; i < termsDual.size(); ++i)
		{
			this->xTilde(data, termsDual[i], dcd[i], dcp[i]);
		}
		//if (doTime) timer.end();
		//if (doTime) printf("Time for termsDual[i]->xTilde(data,tau , dcd[i], dcp[i]); was: %f\n", timer.elapsed());

		//if (doTime) timer.reset();
		for (int i = 0; i < termsPrimal.size(); ++i)
		{
			termsPrimal[i]->applyProx(data, tau, pcp[i]);
		}
		//if (doTime) timer.end();
		//if (doTime) printf("Time for termsPrimal[i]->applyProx(data, tau, pcp[i]); was: %f\n", timer.elapsed());

		//if (doTime) timer.reset();
		//do overrelaxation
		for (int i = 0; i < data->xTilde.size(); ++i)
		{
			doOverrelaxation(data->x[i], data->xOld[i], data->xBar[i]);
		}
		//if (doTime) timer.end();
		//if (doTime) printf("Time for doOverrelaxation(data->x[i], data->xOld[i], data->xBar[i]); was: %f\n", timer.elapsed());


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
		for (int i = 0; i < data->getNumPrimalVars(); ++i)
		{
			calculateXYError(data->x[i], data->xOld[i], data->xError[i], tau[i]);
		}

		for (int i = 0; i < data->getNumDualVars(); ++i)
		{
			calculateXYError(data->y[i], data->yOld[i], data->yError[i], sigma[i]);
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