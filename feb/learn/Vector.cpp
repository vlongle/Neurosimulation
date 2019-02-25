#include "Vector.h"


template <typename T>
void Vector<T>::put(const int index, T elt){
            sum = sum - (*this)[index] + elt;
            (*this)[index] = elt;

}


