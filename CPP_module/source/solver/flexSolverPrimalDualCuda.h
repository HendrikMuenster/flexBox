#ifndef flexSolverPrimalDualCuda_H
#define flexSolverPrimalDualCuda_H

#include "flexSolver.h"

using namespace thrust::placeholders;

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
    
    //List of primal variables corresponding to primal terms
	std::vector<std::vector<int> > pcp;
	//List of primal variables corresponding to dual terms
	std::vector<std::vector<int> > dcp;
	//List of dual variables corresponding to dual terms
	std::vector<std::vector<int> > dcd;

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

		this->calculateTauSigma(data);
	}
	
	void calculateTauSigma(flexBoxData<T, Tvector> *data)
	{
		for (int k = 0; k < (int)termsDual.size(); ++k)
		{
			auto numPrimals = (int)dcp[k].size();
			auto numDuals = (int)dcd[k].size();

			for (int i = 0; i < numDuals; ++i)
			{
				const int dualNum = dcd[k][i];

				for (int j = 0; j < numPrimals; ++j)
				{
					const int primalNum = dcp[k][j];

					auto operatorNumber = numPrimals * i + j;

					auto tmpVec = termsDual[k]->operatorList[operatorNumber]->getAbsRowSumCUDA();
                    thrust::transform(tmpVec.begin(), tmpVec.end(), data->sigmaElt[dualNum].begin(), data->sigmaElt[dualNum].begin(), thrust::plus<T>());

					tmpVec = termsDual[k]->operatorListT[operatorNumber]->getAbsRowSumCUDA();
                    thrust::transform(tmpVec.begin(), tmpVec.end(), data->tauElt[primalNum].begin(), data->tauElt[primalNum].begin(), thrust::plus<T>());
				}
			}
		}

		for (int i = 0; i < (int)data->y.size(); ++i)
		{
            thrust::transform(data->sigmaElt[i].begin(), data->sigmaElt[i].end(), data->sigmaElt[i].begin(), (T)1 / (0.00001 + _1) );
		}

		for (int i = 0; i < data->x.size(); ++i)
		{
            thrust::transform(data->tauElt[i].begin(), data->tauElt[i].end(), data->tauElt[i].begin(), (T)1 / (0.00001 + _1) );
		}
	}

	void addPrimal(flexTermPrimal<T, Tvector>* _primalPart, std::vector<int> _correspondingPrimals)
	{
		termsPrimal.push_back(_primalPart); 
		
		pcp.push_back(_correspondingPrimals);
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

		dcp.push_back(_correspondingPrimals);
		dcd.push_back(tmpDCD);
	}
    
    struct yTildeFunctor
    {
        __host__ __device__
        yTildeFunctor(){};

        template <typename Tuple>
        __host__ __device__
        void operator()(Tuple t)
        {
            thrust::get<0>(t) = thrust::get<1>(t) + thrust::get<2>(t) * thrust::get<0>(t);
        }
    };

	void yTilde(flexBoxData<T, Tvector>* data, flexTermDual<T, Tvector>* dualTerm, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		for (int i = 0; i < (int)dualNumbers.size(); ++i)
		{
			const int dualNum = dualNumbers[i];
            
            //yTilde = 0
            vectorScalarSet(data->yTilde[dualNum],(T)0);
            
			for (int j = 0; j < (int)primalNumbers.size(); ++j)
			{
				int operatorNumber = (int)primalNumbers.size() * i + j;
                
                //yTilde = yTilde + KxBar
				dualTerm->operatorList[operatorNumber]->timesPlus(data->xBar[primalNumbers[j]], data->yTilde[dualNum]);
			}
            
            //yTilde = yOld + sigma * yTilde
            thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(data->yTilde[dualNum].begin(), data->yOld[dualNum].begin(), data->sigmaElt[dualNum].begin())),
                thrust::make_zip_iterator(thrust::make_tuple(data->yTilde[dualNum].end(),   data->yOld[dualNum].end(),   data->sigmaElt[dualNum].end())),
			yTildeFunctor());
		}
	}

	void xTilde(flexBoxData<T, Tvector>* data, flexTermDual<T, Tvector>* dualTerm, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		for (int i = 0; i < (int)dualNumbers.size(); ++i)
		{
			for (int j = 0; j < (int)primalNumbers.size(); ++j)
			{
				int operatorNumber = (int)primalNumbers.size() * i + j;

				const int primalNum = primalNumbers[j];
				const int dualNum = dualNumbers[i];

				dualTerm->operatorListT[operatorNumber]->timesPlus(data->y[dualNum], data->xTilde[primalNum]);
			}
		}
	}
    
    struct xTildeFunctor
    {
        __host__ __device__
        xTildeFunctor(){};

        template <typename Tuple>
        __host__ __device__
        void operator()(Tuple t)
        {
            thrust::get<0>(t) = thrust::get<1>(t) - thrust::get<2>(t) * thrust::get<0>(t);
        }
    };

	void doIteration(flexBoxData<T, Tvector> *data)
	{
        //Timer timer;
        
        //timer.reset();
        data->xOld.swap(data->x);
        data->yOld.swap(data->y);
        //timer.end();printf("Time for swap was: %f\n", timer.elapsed());
        
        //if (doTime) timer.reset();
		for (int i = 0; i < (int)termsDual.size(); ++i)
		{
			//timer.reset();
			this->yTilde(data, termsDual[i], dcd[i], dcp[i]);
			//timer.end(); printf("Time for termsDual[i]->yTilde(data,sigma, dcd[i], dcp[i]); was: %f\n", timer.elapsed());
            
			//timer.reset();
			termsDual[i]->applyProx(data, dcd[i], dcp[i]);
			//timer.end();printf("Time for termsDual[i]->applyProx(data,sigma, dcd[i], dcp[i]); was: %f\n", timer.elapsed());
		}

		//timer.reset();
		for (int i = 0; i < (int)data->xTilde.size(); ++i)
		{
			vectorScalarSet(data->xTilde[i], (T)0);
		}
		//timer.end(); printf("Time for vectorScalarSet(data->xTilde[i], (T)0); was: %f\n", timer.elapsed());

		//timer.reset();
		for (int i = 0; i < (int)termsDual.size(); ++i)
		{
			this->xTilde(data, termsDual[i], dcd[i], dcp[i]);
		}
        //timer.end(); printf("Time for this->xTilde was: %f\n", timer.elapsed());
        
        //timer.reset();
		// set tildeX = xOld - tau * tildeX
		for (int primalNum = 0; primalNum < data->xTilde.size(); ++primalNum)
		{
            thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(data->xTilde[primalNum].begin(), data->xOld[primalNum].begin(), data->tauElt[primalNum].begin())),
                thrust::make_zip_iterator(thrust::make_tuple(data->xTilde[primalNum].end(),   data->xOld[primalNum].end(),   data->tauElt[primalNum].end())),
			xTildeFunctor());
		}
        //timer.end(); printf("Time for this->xTilde2 was: %f\n", timer.elapsed());
        

		//timer.reset();
		for (int primalNum = 0; primalNum < (int)termsPrimal.size(); ++primalNum)
		{
			termsPrimal[primalNum]->applyProx(data, pcp[primalNum]);
		}
		//timer.end(); printf("Time for termsPrimal[i]->applyProx(data, tau, pcp[i]); was: %f\n", timer.elapsed());

		//do overrelaxation
		//timer.reset();
		for (int primalNum = 0; primalNum < (int)data->xTilde.size(); ++primalNum)
		{
            thrust::transform(data->x[primalNum].begin(), data->x[primalNum].end(),data->xOld[primalNum].begin(), data->xBar[primalNum].begin(), _1 + _1 - _2);
		}
        //timer.end(); printf("Time for overrelax was: %f\n", timer.elapsed());
        
        //printf("\n\n\nNewIt\n\n\n");
	}
    
	void yError(flexBoxData<T, Tvector>* data, flexTermDual<T, Tvector>* dualTerm, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		for (int i = 0; i < (int)dualNumbers.size(); ++i)
		{
			for (int j = 0; j < (int)primalNumbers.size(); ++j)
			{
				//Check here
				int operatorNumber = (int)primalNumbers.size() * i + j;

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
		for (int i = 0; i < (int)dualNumbers.size(); ++i)
		{
			for (int j = 0; j < (int)primalNumbers.size(); ++j)
			{
				//Check here
				int operatorNumber = (int)primalNumbers.size() * i + j;

				const int primalNum = primalNumbers[j];
				const int dualNum = dualNumbers[i];

				data->yTmp[dualNum] = data->y[dualNum];
				vectorMinus(data->yTmp[dualNum], data->yOld[dualNum]);

				dualTerm->operatorListT[operatorNumber]->timesMinus(data->yTmp[dualNum], data->xError[primalNum]);
			}
		}
	}
    
    struct calculateErrorFunctor
    {
        __host__ __device__
        calculateErrorFunctor(){};

        template <typename Tuple>
        __host__ __device__
        void operator()(Tuple t)
        {
            thrust::get<0>(t) = (thrust::get<1>(t) - thrust::get<2>(t)) / thrust::get<3>(t);
        }
    };

	T calculateError(flexBoxData<T, Tvector> *data)
	{
		//first part of primal and dual residuals
        //xError = (x - xOld) / tau
		for (int i = 0; i < (int)data->x.size(); ++i)
		{
            thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(data->xError[i].begin(), data->x[i].begin(), data->xOld[i].begin(), data->tauElt[i].begin())),
                thrust::make_zip_iterator(thrust::make_tuple(data->xError[i].end(),   data->x[i].end(),   data->xOld[i].end(),   data->tauElt[i].end())),
            calculateErrorFunctor());
		}

		for (int i = 0; i < (int)data->y.size(); ++i)
		{
            thrust::for_each(
                thrust::make_zip_iterator(thrust::make_tuple(data->yError[i].begin(), data->y[i].begin(), data->yOld[i].begin(), data->sigmaElt[i].begin())),
                thrust::make_zip_iterator(thrust::make_tuple(data->yError[i].end(),   data->y[i].end(),   data->yOld[i].end(),   data->sigmaElt[i].end())),
            calculateErrorFunctor());
		}

		//operator specifiy error
		for (int i = 0; i < (int)termsDual.size(); ++i)
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