/* Clock.cpp
 * Data structures for our simulation
 * Vector class inherits from <vector> and add these methods
 *          put(index, element): put element at vector[index] and update sum
 *          remove_elt(index): more efficient than vector::erase
 *          uniformly_pop(random_engine): select a vector[index] and pop back the element
 */

#include <random>
#include <vector>
#include<iostream>
using namespace std;

#ifndef Vector_H_
#define Vector_H_






template <class T>
class Vector: public vector<T>{
    private:
        T sum;
    public:
    
    void put(const int index, T elt){
        sum = sum - (*this)[index] + elt;
        (*this)[index] = elt;
    }
    
    bool remove_elt(unsigned index){
                  /* return true/ false: success/ failure. */
                    // bound checking
                    if (this->empty() || index < 0 || index >= this -> size()){
                        cout << (this) << " is empty vector" << endl;
                        return false;
                    }
        
                    else{
                        // "move" the last element to the indexed element
                        // effectively erasing the indexed element and shrink the vector by pop_back
                        T elt = (*this)[(*this).size()-1];
                        (*this)[index] = elt;
                        (*this).pop_back();
                        return true;
                    }
        
                }
    T get_sum(){
        return sum;
    }
    
    
    T pop(int index){
                // bound-checking
                if ((*this).size() == 0){
                    cout << "Vector is empty" << endl;
                    return -1;
                }
                else{
        
                    T result = (*this)[index];
                    remove_elt(index);
                    return result;
        
                }
        
    }
//    T uniformly_select(RandomEngine *randomEng){
//
//        // bound-checking
//        if ((*this).size() == 0){
//            cout << "Vector is empty" << endl;
//            return -1;
//        }
//        else{
//
//            int random_index = floor(u(mt)*(*this).size());
//            T result = (*this)[random_index];
//            remove_elt(random_index);
//            return result;
//
//        }
//    }
    
    
    void update_sum() {
        T tmp_sum = 0;
        for (auto&& it : (*this)) tmp_sum += it;
        //       cout<<tmp_sum<<endl;
        sum = tmp_sum;
    }
    
    int select(mt19937& mt, uniform_real_distribution<double>& u)
    {
        if( this->size() == 0)
        {
            cout<<"Size cannot be zero "<<endl;
            return -1;
        }
        else
        {
            int index = floor(u(mt)*this->size());
            T tmp =  (*this)[index];
            remove_elt(index);
            return tmp;
        }
    }

};
#endif
