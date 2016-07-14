#ifndef flexBoxDataGPU_H
#define flexBoxDataGPU_H

#include "flexBoxData.h"

#include <iostream>
#include <vector>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>


template < class T, class Tdata >
class flexBoxDataGPU : public flexBoxData<T, Tdata>
{
public:
	flexBoxDataGPU() : flexBoxData<T, Tdata>(){}

	void addPrimalVar(int numberOfElements)
	{
		Tdata emptyX(numberOfElements, static_cast<T>(0));

		this->x.push_back(emptyX);
		this->xTmp.push_back(emptyX);
		this->xOld.push_back(emptyX);
		this->xBar.push_back(emptyX);
		this->xTilde.push_back(emptyX);
		this->xError.push_back(emptyX);

		this->tauElt.push_back(emptyX);
	}

	void addDualVar(int numberOfElements)
	{
		Tdata emptyY(numberOfElements, static_cast<T>(0));

		this->y.push_back(emptyY);
		this->yTmp.push_back(emptyY);
		this->yOld.push_back(emptyY);
		this->yTilde.push_back(emptyY);
		this->yError.push_back(emptyY);

		this->sigmaElt.push_back(emptyY);
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
		std::vector<T> xTmp(this->x[i].size());

		thrust::copy(this->x[i].begin(), this->x[i].end(), xTmp.begin());

		return xTmp;
	}

	std::vector<T> getDual(int i)
	{
		std::vector<T> yTmp(this->y[i].size());

		thrust::copy(this->y[i].begin(), this->y[i].end(), yTmp.begin());

		return yTmp;
	}

	void setPrimal(int i, std::vector<T> input)
	{
		thrust::copy(input.begin(), input.end(), this->x[i].begin());
		thrust::copy(input.begin(), input.end(), this->xOld[i].begin());
		thrust::copy(input.begin(), input.end(), this->xBar[i].begin());
	}

	void setDual(int i, std::vector<T> input)
	{
		thrust::copy(input.begin(), input.end(), this->y[i].begin());
		thrust::copy(input.begin(), input.end(), this->yOld[i].begin());
		thrust::copy(input.begin(), input.end(), this->yTilde[i].begin());
	}
};

#endif