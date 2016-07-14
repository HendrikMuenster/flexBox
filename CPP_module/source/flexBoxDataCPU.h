#ifndef flexBoxDataCPU_H
#define flexBoxDataCPU_H

#include "flexBoxData.h"

template < class T, class Tdata >
class flexBoxDataCPU : public flexBoxData<T, Tdata>
{
	public:

	flexBoxDataCPU() : flexBoxData<T, Tdata>(){}

	void addPrimalVar(int numberOfElements)
	{
		Tdata emptyX(numberOfElements, static_cast<T>(0));

		this->x.push_back(emptyX);
		this->xTmp.push_back(emptyX);
		this->xOld.push_back(emptyX);
		this->xBar.push_back(emptyX);
		this->xTilde.push_back(emptyX);
		this->xError.push_back(emptyX);
	}

	void addDualVar(int numberOfElements)
	{
		Tdata emptyY(numberOfElements, static_cast<T>(0));

		this->y.push_back(emptyY);
		this->yTmp.push_back(emptyY);
		this->yOld.push_back(emptyY);
		this->yTilde.push_back(emptyY);
		this->yError.push_back(emptyY);
	}

	int getNumPrimalVars()
	{
		return this->x.size();
	}

	int getNumDualVars()
	{
		return this->y.size();
	}

	std::vector<T> getPrimal(int i)
	{
		return this->x[i];
	}

	std::vector<T> getDual(int i)
	{
		return this->y[i];
	}

	void setPrimal(int i, std::vector<T> input)
	{
		this->x[i] = input;
	}

	void setDual(int i, std::vector<T> input)
	{
		this->y[i] = input;
	}
};

#endif