#ifndef RandomEngine_H_
#define RandomEngine_H_

#include <random>
#include "Vector.h"
using namespace std;

// Each population should have one separate RandomEngine.
struct RandomGenerator{
    mt19937 mt;
    uniform_real_distribution<double> u;
};

class RandomEngine{
//private:
//    struct RandomGenerator randGen;
public:
    struct RandomGenerator randGen;
    void uniform_in_range(int range);
    void set_random_generator(mt19937 set_mt, uniform_real_distribution<double> set_u);
    
    
    
    // return index in an array based on the "weight" of element stored at each location.
    int weighted_random_select(Vector<double>* array); // used to randomly choose an action to take during each timestep.
    
    double get_random_number();
    
};


#endif
