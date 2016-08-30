


#ifndef flexbox_H
#define flexbox_H

#include "flexTermDual.h"
#include "flexDualizedOperator.h"
#include "flexDualizedDataTerm.h"

#include "flexTermPrimal.h"

#include "flexSolver.h"

#if __CUDACC__
	#include "flexSolverPrimalDualCuda.h"
	#include "flexBoxDataGPU.h"

	#include <device_functions.h>
#else
	#include "flexSolverPrimalDual.h"
#endif

#include "flexBoxDataCPU.h"



#include <vector>

template < typename T, typename Tdata >
class flexBox
{
	typedef flexLinearOperator < T, Tdata > linOpType;

	private:

	public:
		T tol;

		//Maximum number of Iterations
		int maxIterations;
		int checkError;
		int displayStatus;
		int verbose;

		//List of dimensions
		std::vector<std::vector<int> > dims;
		
		flexBoxData<T,Tdata>* data;
		flexSolver<T, Tdata>* solver;

		bool isMATLAB; // indicates whether flexBox is used via MALTAB

		flexBox(void)
		{
			this->tol = static_cast<T>(1e-5);

			this->maxIterations = static_cast<int>(10000);
			this->checkError = static_cast<int>(100);
			this->displayStatus = static_cast<int>(1000);
			this->verbose = static_cast<int>(0);

			#if __CUDACC__
				this->data = new flexBoxDataGPU<T, Tdata>();
				this->solver = new flexSolverPrimalDualCuda<T, Tdata>();
			#else
				this->data = new flexBoxDataCPU<T, Tdata>();
				this->solver = new flexSolverPrimalDual<T, Tdata>();
			#endif	

			this->isMATLAB = false;

		};

		~flexBox()
		{
			delete data;
			delete solver;
		}

		int getNumPrimalVars() const 
		{
			return data->getNumPrimalVars();
		}

		int getNumDualVars() const
		{
			return data->getNumDualVars();
		}

		std::vector<T> getPrimal(int i)
		{
			return data->getPrimal(i);
		}

		std::vector<T> getDual(int i)
		{
			return data->getDual(i);
		}

		void setPrimal(int i, std::vector<T> input)
		{
			data->setPrimal(i, input);
		}

		void setDual(int i, std::vector<T> input)
		{
			data->setDual(i, input);
		}

		void init()
		{
			solver->init();
		}

		std::vector<int> getDims(int i)
		{
			return dims.at(i);
		}

		int addPrimalVar(std::vector<int> _dims)
		{
			int numberOfElements = vectorProduct(_dims);

			data->addPrimalVar(vectorProduct(_dims));

			dims.push_back(_dims);

			return getNumPrimalVars() - 1;
		}

		void addPrimal(flexTermPrimal<T,Tdata>* _primalPart, std::vector<int> _correspondingPrimals)
		{
			solver->addPrimal(_primalPart, _correspondingPrimals);
		}

		void addDual(flexTermDual<T, Tdata>* _dualPart, std::vector<int> _correspondingPrimals)
		{
			solver->addDual(data, _dualPart, _correspondingPrimals);
		}

		void runAlgorithm()
		{
			solver->init(data);

			T error = static_cast<int>(1);
			int iteration = 0;
			
			Timer timer;
			Timer timer2;
			bool doTime = true;

			if (doTime) timer.reset();

			while (error > tol && iteration < maxIterations)
			{
				//timer2.reset();
				solver->doIteration(data);
				//timer2.end(); printf("Time for iteration was: %f\n", timer2.elapsed());
				
				if (iteration % displayStatus == 1)
				{
					if (this->isMATLAB)
					{
						if (this->verbose > 0)
						{
							mexPrintf("Iteration #%d | Error:%f\n", iteration, error);
							mexEvalString("pause(.0001);");
						}
					}
					else
					{

					}
				}
				
				if (iteration % checkError == 0)
				{
					error = solver->calculateError(data);
				}

				++iteration;
			}

			if (doTime) timer.end();
			
			if (this->verbose > 0)
			{
				if (doTime) printf("Time for %d Iterations was: %f\n", iteration, timer.elapsed());
			}

		}
};

#endif