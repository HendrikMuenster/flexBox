#ifndef flexProx_H
#define flexProx_H



#include <vector>

template < class T, class Tvector >
class flexProx
{
private:

public:

	flexProx()
	{
		
	}

	~flexProx()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}

	virtual void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers) = 0;
	
	virtual void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<std::vector<T>> &fList) = 0;
};

#endif