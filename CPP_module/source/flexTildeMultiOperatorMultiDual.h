#ifndef flexTildeMultiOperatorMultiDual_H
#define flexTildeMultiOperatorMultiDual_H

#include "flexBox.h"

template < class T, class Tvector >
class flexTildeMultiOperatorMultiDual : public flexTermDual<T, Tvector>
{
private:


public:
	flexTildeMultiOperatorMultiDual(prox _p,T _alpha, int numberPrimals, int _numberVars) : flexTermDual<T,Tvector>(_p,_alpha, numberPrimals, _numberVars){};


};

#endif