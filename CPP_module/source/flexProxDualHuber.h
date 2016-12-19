#ifndef flexProxDualHuber_H
#define flexProxDualHuber_H



#include "flexProx.h"

template < typename T, typename Tvector >
class flexProxDualHuber : public flexProx<T, Tvector>
{
private:
	T huberEpsilon;
public:

	flexProxDualHuber(T _huberEpsilon) : flexProx<T, Tvector>()
	{
		huberEpsilon = _huberEpsilon;
	}

	~flexProxDualHuber()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}

	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers)
	{
		if (dualNumbers.size() == 1)
		{
			T* ptrY0 = data->y[dualNumbers[0]].data();

			T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();

			T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

			int numElements = data->yTilde[dualNumbers[0]].size();

			#pragma omp parallel for
			for (int i = 0; i < numElements; i++)
			{
				T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * (T)this->huberEpsilon / alpha);

				T yTmp = (T)1 / myMax<T>((T)1, fabs(ptrYtilde0[i] * huberFactor) / alpha);

				ptrY0[i] = ptrYtilde0[i] * yTmp;
			}
		}
		else if (dualNumbers.size() == 2)
		{
			T* ptrY0 = data->y[dualNumbers[0]].data();
			T* ptrY1 = data->y[dualNumbers[1]].data();

			T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
			T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();

			T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

			int numElements = data->yTilde[dualNumbers[0]].size();

			#pragma omp parallel for
			for (int i = 0; i < numElements; i++)
			{
				T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * (T)this->huberEpsilon / alpha);

				T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i] * huberFactor) + pow2(ptrYtilde1[i] * huberFactor)) / alpha);

				ptrY0[i] = ptrYtilde0[i] * yTmp;
				ptrY1[i] = ptrYtilde1[i] * yTmp;
			}
		}
		else if (dualNumbers.size() == 3)
		{
			T* ptrY0 = data->y[dualNumbers[0]].data();
			T* ptrY1 = data->y[dualNumbers[1]].data();
			T* ptrY2 = data->y[dualNumbers[2]].data();

			T* ptrYtilde0 = data->yTilde[dualNumbers[0]].data();
			T* ptrYtilde1 = data->yTilde[dualNumbers[1]].data();
			T* ptrYtilde2 = data->yTilde[dualNumbers[2]].data();

			T* ptrSigma = data->sigmaElt[dualNumbers[0]].data();

			int numElements = data->yTilde[dualNumbers[0]].size();

			#pragma omp parallel for
			for (int i = 0; i < numElements; i++)
			{
				T huberFactor = (T)1 / ((T)1.0 + ptrSigma[i] * (T)this->huberEpsilon / alpha);

				T yTmp = (T)1 / myMax<T>((T)1, sqrtf(pow2(ptrYtilde0[i] * huberFactor) + pow2(ptrYtilde1[i] * huberFactor) + pow2(ptrYtilde2[i] * huberFactor)) / alpha);

				ptrY0[i] = ptrYtilde0[i] * yTmp;
				ptrY1[i] = ptrYtilde1[i] * yTmp;
				ptrY2[i] = ptrYtilde2[i] * yTmp;
			}
		}
		else
		{
			printf("Alert! Huber prox not implemented for dim>3");
		}
	}
	
	void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<std::vector<T>> &fList)
	{

	}
};

#endif