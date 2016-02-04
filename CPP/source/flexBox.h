#ifndef flexbox_H
#define flexbox_H

#include "flexVector.h"
#include "flexTermDual.h"
#include "flexTermPrimal.h"
#include "flexBoxData.h"



template < typename T >
class flexBox
{
	private:
		//Params
		T tol;
		T theta;
		flexVector<T> tau;
		flexVector<T> sigma;

		//Maximum number of Iterations
		int maxIterations;
		int numberPrimalVars;
		int checkError;

		//List of primal terms is list of pointers to terms
		flexVector<flexTermPrimal<T>*> termsPrimal;
		//List of dual terms is list of pointers to terms
		flexVector<flexTermDual<T>*> termsDual;

		flexBoxData<T> data;

		//List of dimensions
		flexVector<flexVector<int>> dims;
		//List of primal variables corresponding to primal terms
		flexVector<flexVector<int>> pcp;
		//List of primal variables corresponding to dual terms
		flexVector<flexVector<int>> dcp;
		//List of dual variables corresponding to dual terms
		flexVector<flexVector<int>> dcd;

	public:

		bool isMATLAB; // indicates whether flexBox is used via MALTAB

		flexBox(void)
		{
			tol = static_cast<T>(1e-5);
			theta = static_cast<T>(1);


			maxIterations = static_cast<int>(1000);
			numberPrimalVars = static_cast<int>(0);
			checkError = static_cast<int>(100);

			isMATLAB = false;
		};

		flexVector<T> getPrimal(int i)
		{
			return data.x[i];
		}

		void init()
		{
			//zero init
			for (int i = 0; i < data.x.size(); ++i)
			{
				tau.push_back(static_cast<T>(0));
			}

			for (int i = 0; i < data.y.size(); ++i)
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
				tau[i] = 1 / tau[i];
			}

			for (int i = 0; i < sigma.size(); ++i)
			{
				sigma[i] = 1 / sigma[i];
			}
		}

		flexVector<int> getDims(int i)
		{
			return dims[i];
		}

		void addPrimal(flexTermPrimal<T>* _primalPart)
		{
			termsPrimal.push_back(_primalPart);

			int numberOfElements = _primalPart->getDims().product();

			flexVector<T> emptyX(numberOfElements,static_cast<T>(0));

			flexVector<int> tmpPCP;
			for (int i = 0; i < _primalPart->getNumberVars(); ++i)
			{
				data.x.push_back(emptyX);
				data.xOld.push_back(emptyX);
				data.xBar.push_back(emptyX);
				data.xTilde.push_back(emptyX);
				data.xError.push_back(emptyX);

				flexVector<int> tmpDims = _primalPart->getDims();
				dims.push_back(tmpDims);

				tmpPCP.push_back(numberPrimalVars);

				++numberPrimalVars;
			}


			pcp.push_back(tmpPCP);
		}

		void addDual(flexTermDual<T>* _dualPart,flexVector<int> _correspondingPrimals)
		{
			termsDual.push_back(_dualPart);

			flexVector<int> tmpDCD;
			for (int i = 0; i < _dualPart->getNumberVars(); ++i)
			{
				flexVector<T> emptyY(_dualPart->dualVarLength(i), static_cast<T>(0));

				tmpDCD.push_back(data.y.size());

				data.y.push_back(emptyY);
				data.yOld.push_back(emptyY);
				data.yTilde.push_back(emptyY);
				data.yError.push_back(emptyY);
			}

			dcp.push_back(_correspondingPrimals);
			dcd.push_back(tmpDCD);
		}

		void doIteration()
		{
			//data.x[0].print();

			for (int i = 0; i < data.xOld.size(); ++i)
			{
				data.xOld[i] = data.x[i];
			}

			for (int i = 0; i < data.yOld.size(); ++i)
			{
				data.yOld[i] = data.y[i];
			}

			for (int i = 0; i < termsDual.size(); ++i)
			{
				termsDual[i]->yTilde(data,sigma, dcd[i], dcp[i]);
				termsDual[i]->applyProx(data,sigma, dcd[i], dcp[i]);
			}

			//data.yTilde[0].print();

			for (int i = 0; i < data.xTilde.size(); ++i)
			{
				data.xTilde[i] = data.x[i];
			}

			for (int i = 0; i < termsDual.size(); ++i)
			{
				termsDual[i]->xTilde(data,tau , dcd[i], dcp[i]);
			}

			for (int i = 0; i < termsPrimal.size(); ++i)
			{
				termsPrimal[i]->applyProx(data, tau, pcp[i]);
			}

			/*data.xBar = data.x;
			data.xBar += data.x;
			data.xBar -= data.xOld;*/

			for (int i = 0; i < data.xTilde.size(); ++i)
			{
				//musst be implemented in the manner below, times theta missing
				//data.xBar[i] = data.x[i] + (data.x[i] - data.xOld[i]);
				data.xBar[i] = data.x[i];
				data.xBar[i] += data.x[i];
				data.xBar[i] -= data.xOld[i];
			}

			//data.xBar[0].print();
			//data.x[0].print();

		}

		T calculateError(void)
		{
			T error = static_cast<int>(1);

			//first part of primal and dual residuals
			for (int i = 0; i < data.xOld.size(); ++i)
			{
				flexVector<T> tmpVec(data.x[i]);
				tmpVec -= data.xOld[i];
				tmpVec /= tau[i];

				data.xError[i] = tmpVec;
			}

			for (int i = 0; i < data.yOld.size(); ++i)
			{
				flexVector<T> tmpVec(data.y[i]);
				tmpVec -= data.yOld[i];
				tmpVec /= sigma[i];

				data.yError[i] = tmpVec;
			}

			//operator specifiy error
			for (int i = 0; i < termsDual.size(); ++i)
			{
				termsDual[i]->xError(data, tau, dcd[i], dcp[i]);
				termsDual[i]->yError(data, sigma, dcd[i], dcp[i]);
			}

			//sum things up
			T primalResidual = static_cast<T>(0);
			T dualResidual = static_cast<T>(0);

			for (int i = 0; i < data.xOld.size(); ++i)
			{
				data.xError[i].abs(); //take absolute value of entries

				primalResidual += data.xError[i].sum() / data.xError[i].size(); //sum values up and add to primal residual
			}
			primalResidual = primalResidual / data.xOld.size();

			for (int i = 0; i < data.yOld.size(); ++i)
			{
				data.yError[i].abs(); //take absolute value of entries

				dualResidual += data.yError[i].sum() / data.yError[i].size(); //sum values up and add to dual residual
			}
			dualResidual = dualResidual / data.yOld.size();

			error = primalResidual + dualResidual;

			return error;
		}

		void runAlgorithm()
		{
			init();

			T error = static_cast<int>(1);
			int iteration = 1;

			while (error > tol && iteration < maxIterations)
			{
				doIteration();
				//printf("Value f at 5 is %f\n",data.x[0][5]);
				if (iteration % 100 == 1)
				{
					if (isMATLAB == true)
					{
						printf("Iteration #%d | Error:%f\n", iteration, error);
						mexEvalString("drawnow;");
					}
					else
					{

					}
					
				}
				
				if (iteration % checkError == 1)
				{
					error = calculateError();
				}
				

				++iteration;
			}

		}
};

#endif