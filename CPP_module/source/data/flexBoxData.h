#ifndef flexBoxData_H
#define flexBoxData_H

#include <vector>

template < class T, class Tdata >
class flexBoxData
{
	public:
	//List of primal variables
	std::vector<Tdata> x, xTmp, xOld, xTilde, xBar, xError;
	//List of dual variables
	std::vector<Tdata> y, yTmp, yOld, yTilde, yError;
	//weights
	std::vector<Tdata> tauElt;
	std::vector<Tdata> sigmaElt;

	flexBoxData(){}


	virtual void addPrimalVar(int numberOfElements) = 0;
	virtual void addDualVar(int numberOfElements) = 0;

	virtual int getNumPrimalVars() = 0;
	virtual int getNumDualVars() = 0;

	virtual std::vector<T> getPrimal(int i) = 0;
	virtual std::vector<T> getDual(int i) = 0;

	virtual void setPrimal(int i, std::vector<T> input) = 0;
	virtual void setDual(int i, std::vector<T> input) = 0;
};

#endif