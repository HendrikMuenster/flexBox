#ifndef flexProxList_H
#define flexProxList_H

#include "flexBoxData.h"
#include "vector"

template < typename T, typename Tvector>
class flexProxList
{

public:

	flexProxList() {};

	static void primalEmptyProx(flexBoxData<T, Tvector>* data, std::vector<T> tau, std::vector<int> primalNumbers)
	{
		for (int i = 0; i < primalNumbers.size(); ++i)
		{
			data->x[primalNumbers[i]].swap(data->xTilde[primalNumbers[i]]);
		}
	};

	static void dualL1AnisoProx(flexBoxData<T, Tvector>* data, const std::vector<T> &sigma, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, T alpha)
	{
		for (int j = 0; j < dualNumbers.size(); ++j)
		{
			vectorProjectL1Aniso<T>(data->yTilde[dualNumbers[j]], data->y[dualNumbers[j]], alpha);
		}
	};

	static void dualL1IsoProx(flexBoxData<T, Tvector>* data, const std::vector<T> &sigma, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, T alpha)
	{
		if (dualNumbers.size() == 2)
		{
			vectorProjectL1Iso2D<T>(data->yTmp[dualNumbers[0]], data->yTilde[dualNumbers[0]], data->yTilde[dualNumbers[1]], data->y[dualNumbers[0]], data->y[dualNumbers[1]], alpha);
		}
		else
		{
			printf("Alert! Iso prox not implemented for dim>2");
		}
	};

	static void dualL2Prox(flexBoxData<T, Tvector>* data, const std::vector<T> &sigma, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, T alpha)
	{
		for (int j = 0; j < dualNumbers.size(); ++j)
		{
			const int dualNum = dualNumbers[j];
			vectorProjectL2<T>(data->yTilde[dualNum], data->y[dualNum], alpha, sigma[dualNum]);
		}
	};

	static void dualFrobeniusProx(flexBoxData<T, Tvector>* data, const std::vector<T> &sigma, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, T alpha)
	{
		if (dualNumbers.size() == 2)
		{
			vectorProjectFrobenius2D<T>(data->yTmp[dualNumbers[0]], data->yTilde[dualNumbers[0]], data->yTilde[dualNumbers[1]], data->y[dualNumbers[0]], data->y[dualNumbers[1]], alpha);
		}
		else
		{
			printf("Alert! Frob prox not implemented for dim>2");
		}
	};
	
	static void dualL2DataProx(flexBoxData<T, Tvector>* data, const std::vector<T> &sigma, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, T alpha, Tvector &f)
	{
		vectorProjectL2data<T>(data->y[dualNumbers[0]], data->yTilde[dualNumbers[0]], f, alpha, sigma[dualNumbers[0]]);
	};

	static void dualL1DataProx(flexBoxData<T, Tvector>* data, const std::vector<T> &sigma, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, T alpha, Tvector &f)
	{
		vectorProjectL1data<T>(data->y[dualNumbers[0]], data->yTilde[dualNumbers[0]], f, alpha, sigma[dualNumbers[0]]);
	};

	static void dualKLDataProx(flexBoxData<T, Tvector>* data, const std::vector<T> &sigma, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, T alpha, Tvector &f)
	{
		vectorProjectKLdata<T>(data->y[dualNumbers[0]], data->yTilde[dualNumbers[0]], f, alpha, sigma[dualNumbers[0]]);
	};
};

#endif