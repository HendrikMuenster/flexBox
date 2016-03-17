#ifndef flexBoxData_H
#define flexBoxData_H

#include "flexVector.h"

template < typename T >
class flexBoxData
{
	public: 
	//List of primal variables
	flexVector<flexVector<T>> x, xOld, xTilde, xBar, xError;
	//List of dual variables
	flexVector<flexVector<T>> y, yOld, yTilde, yError;

	flexBoxData(void)
	{

	}
};

#endif