#ifndef flexProx_H
#define flexProx_H


#include "tools.h"
#include <vector>

template < class T, class Tvector >
class flexProx
{
private:
    
public:
    const prox p;
    
	flexProx(prox _p) : p(_p)
	{
		
	}

	~flexProx()
	{
		if (VERBOSE > 0) printf("Destructor prox\n!");
	}
    
    prox getProx()
    {
        return p;
    }

	virtual void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers) = 0;
	
	virtual void applyProx(T alpha, flexBoxData<T, Tvector>* data, const std::vector<int> &dualNumbers, const std::vector<int> &primalNumbers, std::vector<std::vector<T>> &fList) = 0;
};

#endif