#ifndef flexSolver_H
#define flexSolver_H



#include <vector>

template < class T, class Tvector >
class flexSolver
{
private:

public:

	flexSolver()
	{
		
	}

	virtual ~flexSolver()
	{
		if (VERBOSE > 0) printf("Destructor solver\n!");
	}

	virtual void init(flexBoxData<T, Tvector> *data) = 0;

	virtual void addPrimal(flexTermPrimal<T, Tvector>* _primalPart, std::vector<int> _correspondingPrimals) = 0;

	virtual void addDual(flexBoxData<T, Tvector> *data, flexTermDual<T, Tvector>* _dualPart, std::vector<int> _correspondingPrimals) = 0;

	virtual void doIteration(flexBoxData<T, Tvector> *data) = 0;

	virtual T calculateError(flexBoxData<T, Tvector> *data) = 0;
};

#endif