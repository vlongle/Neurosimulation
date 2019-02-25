#include <vector>
#ifndef Vector_H_
#define Vector_H_

template <class T>
class Vector: public std:: vector<T>{
    private: 
        T sum;
    public:
        void put(const int index, T elt);
};



#endif
